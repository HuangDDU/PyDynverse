# TODO: 后续通过继承关系AnnData实现Wrapper
# TODO: 参考dyncli中关于read与write中装饰器的使用方式, 改善代码结构
import os
import json
import subprocess

import pandas as pd
import rpy2.robjects as ro

from .._logging import logger


class DynverseDockerInput():
    def __init__(self, expression, expression_id, cell_ids, feature_ids, parameters, priors, seed, verbose):
        self.expression = expression
        self.expression_id = expression_id
        self.cell_ids = cell_ids
        self.feature_ids = feature_ids
        self.parameters = parameters
        self.priors = priors
        self.seed = seed
        self.verbose = verbose

    def save_json(self, input_json_filename):
        self.input_json_filename = input_json_filename
        # 稀疏矩阵不能JSON序列化，需要手动提取行列索引
        expression_dict = {
            "x": self.expression.data.tolist(),
            "i": self.expression.indices.tolist(),
            "p": self.expression.indptr.tolist(),
            "Dim": self.expression.shape,
            # TODO: 表达矩阵行、列名称，暂时这样写，后续再把真实的拉进来
            "rownames": list(self.cell_ids),
            "colnames": list(self.feature_ids)
        }
        input_json = {
            self.expression_id: expression_dict,
            "expression_id": self.expression_id,
            "parameters": self.parameters,
            "priors": self.priors,
            "seed": self.seed,
            "verbose": self.verbose
        }
        # logger.debug(input_json)
        # 对于稀疏矩阵的特殊处理
        with open(input_json_filename, "w") as f:
            json.dump(input_json, f)
        logger.debug(f"Save json successfully, path: {input_json_filename}")

    def json2h5(self, input_h5_filename):
        self.input_h5_filename = input_h5_filename
        # 调用R脚本，把生成的json文件转化为h5文件，作为dynverse docker容器需要的的输入
        Rscript_filename = f"{os.path.dirname(__file__)}/../rscript/docker_input_json2h5.R"
        logger.debug(f"h52json script: {Rscript_filename}")
        command_list = [Rscript_filename, "--input_json_filename",
                        self.input_json_filename, "--input_h5_filename", input_h5_filename]
        result = subprocess.run(command_list, capture_output=True, text=True)
        logger.debug(result)
        if result.returncode == 0:
            logger.debug("json2h5 successful!")
        else:
            logger.debug("json2h5 failed!")

    def __str__(self) -> str:
        return f"{self.expression}"  # 目前查看稀疏矩阵是最直观的输入


class DynverseDockerOutput():
    def __init__(self):
        self.id = None
        self.pseudotime = None

    def h52json(self, output_h5_filename, output_json_filename):

        self.output_h5_filename = output_h5_filename
        self.output_json_filename = output_json_filename
        # 调用R脚本，把dynverse docker的输出的h5文件转化JSON文件
        Rscript_filename = f"{os.path.dirname(__file__)}/../rscript/docker_output_h52json.R"
        logger.debug(f"h52json script: {Rscript_filename}")
        command_list = [Rscript_filename, "--output_h5_filename",
                        output_h5_filename, "--output_json_filename", output_json_filename]
        result = subprocess.run(command_list, capture_output=True, text=True)
        logger.debug(result)
        if result.returncode == 0:
            logger.debug("h52json successful!")
        else:
            logger.debug("h52json failed!")

    def load_json(self):
        # 读取json
        with open(self.output_json_filename, "r") as f:
            output_json = json.load(f)
        logger.debug(
            f"Save json successfully, path: {self.output_json_filename}")
        # 解析JSON
        # 暂时简单设置属性, 后续需要使用数据类型, 再对应转换
        self.id = output_json["id"]
        self.cell_ids = output_json["cell_ids"]
        # 自动提取JSON中的键值对
        for k, v in output_json.items():
            self.__setattr__(k, v)
        # 对共有的部分属性的数据结构修改, 方便后续调用, 这部分其实就是wrapper对于轨迹推断输出结果的封装
        self.milestone_percentages = pd.DataFrame(self.milestone_percentages)
        self.progressions = pd.DataFrame(self.progressions)
        self.dimred = pd.DataFrame(self.dimred, index=self.cell_ids)
        self.dimred_segment_progressions = pd.DataFrame(
            output_json["dimred_segment_progressions"])
        self.dimred_segment_points = pd.DataFrame(
            output_json["dimred_segment_points"])
        
        # json文件添加，方便查可能不同轨迹推断类型对于wrapper的输出
        self.output_json = output_json
        
    def __str__(self) -> str:
        return f"""
            id: {self.id}, 
            trajectory_type: {self.trajectory_type}, 
            attribute_list: {self.__dict__.keys()}
            """

    def __getitem__(self, key):
        # 通过键名访问属性
        if hasattr(self, key):
            return getattr(self, key)
        else:
            raise KeyError(
                f"'{self.__class__.__name__}' object has no attribute '{key}'")

    def __setitem__(self, key, value):
        # 通过键名设置属性
        setattr(self, key, value)

    def __contains__(self, item):
        return hasattr(self, item)


def write_h5(x, h5_filename, via_json=True):
    expression_id = x["expression_id"]
    if via_json:
        input_json_filename = f"{h5_filename[:-3]}.json"  # 中间json文件
        input_h5_filename = h5_filename
        dynverse_docker_input = DynverseDockerInput(
            expression=x[expression_id], # 从AnnData里提取
            expression_id = expression_id,
            cell_ids=x["cell_ids"],
            feature_ids=x["feature_ids"],
            parameters=x["parameters"],
            priors=x["priors"],
            seed=x["seed"],
            verbose=x["verbose"]
        )
        dynverse_docker_input.save_json(input_json_filename)
        dynverse_docker_input.json2h5(input_h5_filename)
    else:
        # TODO: 直接通过装饰器的自动转换，不使用json交换文件和单独的R脚本
        task = x
        task[expression_id] = None
        ro.globalenv["task"] = ro.ListVector(task) # 待添加的内容转换到R变量里
        ro.globalenv["h5_filename"] = ro.ListVector(h5_filename) # 待添加的内容转换到R变量里
        ro.r("dynutils::write_h5(task, file.path(paths$dir_dynwrap, h5_filename))") # 调用R修改
    


def read_h5(h5_filename, via_json=True):
    if via_json:
        output_h5_filename = h5_filename
        output_json_filename = f"{h5_filename[:-3]}.json"  # 中间json文件
        dynverse_docker_output = DynverseDockerOutput()
        dynverse_docker_output.h52json(output_h5_filename, output_json_filename)
        dynverse_docker_output.load_json()
        return dynverse_docker_output
    else:
        # TODO: 直接通过装饰器的自动转换，不使用json交换文件和单独的R脚本
        return 
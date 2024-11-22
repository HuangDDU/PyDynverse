import time
import tempfile

from .._logging import logger
from .method_extract_args import _method_extract_inputs, _method_extract_priors
from .method_process_definition import get_default_parameters
from .method_create_ti_method_container import \
    _method_execution_preproc_container, \
    _method_execution_execute_container, \
    _method_execution_postproc_container


def _method_execute(
        dataset,  # AnnData
        method,
        parameters,
        give_priors,
        seed,
        verbose,
        return_verbose,
        debug
):
    if debug:
        verbose = True

    timings = {"execution_start": time.time()}

    # 提取模型输入input，这里暂时直接返回AnnData
    inputs = _method_extract_inputs(dataset, method)
    priors = _method_extract_priors(dataset, method)  # 提取先验知识prior
    parameters = get_default_parameters(method)  # 提取函数的默认参数

    # 丢弃关于输入输出流的额外处理

    # 创建临时工作目录，/tmp目录下
    with tempfile.TemporaryDirectory() as tmp_wd:
        logger.debug(f"Temp wd: {tmp_wd}")

        # 容器执行
        timings["method_beforepreproc"] = time.time()
        preproc_meta = _method_execution_preproc_container(
            method,
            inputs,
            tmp_wd,
            priors,
            parameters,
            verbose,
            seed,
            debug
        )  # 构造DynverseDockerInput对象
        timings["method_afterpostproc"] = time.time()

        trajectory = _method_execution_execute_container(
            method=method,
            preproc_meta=preproc_meta,
            tmp_wd=tmp_wd
        )
        timings["method_afterpostproc"] < - time.time()
        _method_execution_postproc_container(preproc_meta=preproc_meta)
        timings["execution_stop"] = time.time()

    output = {
        "summary": "",
        "trajectory": trajectory
    }

    return output

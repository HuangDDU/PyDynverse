from .wrap_add_expression import get_expression

def _method_extract_inputs(dataset, inputs):
    # 提取模型输入input, 表达矩阵
    input_ids = inputs["input_id"][inputs["type"]=="expression"]
    inputs = {}
    for expression_id in input_ids:
        inputs[expression_id] = get_expression(dataset, expression_id)
    inputs["expression_id"] = input_ids[0] # 作为主要的表达矩阵, 例如Component1和Slingshot需要expression, 而monocle_ddrtree需要counts
    # 额外添加细胞名和特证名
    inputs["cell_ids"] =  dataset["cell_ids"]
    inputs["feature_ids"] =  dataset["feature_ids"]
    return inputs

def _method_extract_priors(dataset, inputs, give_prioirs=None):
    # TODO: 提取先验知识prior
    optional_prior_ids = inputs["input_id"][inputs["type"]=="prior_information"]
    priors = None
    return priors
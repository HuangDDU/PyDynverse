from .wrap_add_expression import get_expression
from .._logging import logger


def _method_extract_inputs(dataset, inputs):
    # 提取模型输入input, 表达矩阵
    input_ids = inputs["input_id"][inputs["type"] == "expression"].tolist()
    inputs = {}
    for expression_id in input_ids:
        inputs[expression_id] = get_expression(dataset, expression_id)
    # 作为主要的表达矩阵, 例如Component1和Slingshot需要expression, 而monocle_ddrtree需要counts
    inputs["expression_id"] = input_ids[0]
    # 额外添加细胞名和特证名
    inputs["cell_ids"] = dataset["cell_ids"]
    inputs["feature_ids"] = dataset["feature_ids"]
    return inputs


def _method_extract_priors(
        dataset,
        inputs,
        give_prioirs=None):
    #
    priors = dataset.get("prior_information", {})
    priors["dataset"] = dataset # 这句话感觉可以不用, 有些多余了
    priors_key_list = priors.keys()
    priors_key_set = set(priors_key_list)

    # 必要的先验输入
    required_prior_ids = inputs["input_id"][inputs["required"]&(inputs["type"]=="prior_information")].tolist()
    required_prior_ids_set = set(required_prior_ids)
    if not (required_prior_ids_set <= priors_key_set):
        # 如果必要的先验输入不完全,则报错
        missing_priors = required_prior_ids_set - priors_key_set
        raise Exception(f"""\n
                ! Prior information {' ,'.join(missing_priors)} is missing from dataset {dataset['id']} but is required by the method. \n
                -> If known, you can add this prior information using add_prior_information(dataset, {' ,'.join([str(i)+' = <prior>' for i in missing_priors])}). \n
                -> Otherwise, this method cannot be used.
                     """)
    required_prior = {k:priors[k] for k in required_prior_ids}

    # 可选的先验输入
    optional_prior_ids = inputs["input_id"][(~inputs["required"])&(inputs["type"]=="prior_information")].tolist()
    optional_prior_ids_set = set(optional_prior_ids)
    if not (optional_prior_ids_set <= priors_key_set):
        # 可选的先验输入不完全, 只警告
        missing_priors = list(optional_prior_ids_set - priors_key_set)
        logger.warning(f"""\n
                       Prior information {','.join(missing_priors)} is optional, but missing from dataset {dataset['id']}. \n
                       Will not give this prior to method.
                       """)
    optional_prior = {k:priors[k] for k in list(optional_prior_ids_set & priors_key_set)}

    # 上述这些操作,一方面是判断先验知识是否给够了,另一方面是过滤掉一些无用的先验知识
    priors = required_prior | optional_prior

    return priors

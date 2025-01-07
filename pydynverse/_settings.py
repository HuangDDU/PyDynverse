# 参考
# Scanpy: https://github.com/scverse/scanpy/blob/main/src/scanpy/_settings.py
# OmicVerse: https://github.com/Starlitnightly/omicverse/blob/master/omicverse/_settings.py
class PyDynverseConfig:
    def __init__(self):
        self.backend = None

    def __getitem__(self, key):
        if hasattr(self, key):
            return getattr(self, key)
        else:
            return None

    def __setitem__(self, key, value):
        setattr(self, key, value)


settings = PyDynverseConfig()  # init导入该配置

# 参考：https://github.com/theislab/moscot/blob/main/src/moscot/_logging.py
import logging

__all__ = ["logger"]

def _setup_logger() -> "logging.Logger":
    from rich.console import Console
    from rich.logging import RichHandler

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    console = Console(force_terminal=True)
    if console.is_jupyter is True:
        console.is_jupyter = False
    ch = RichHandler(show_path=False, console=console)
    # ch.setFormatter(logging.Formatter("%(name)s -%(filename)s:%(lineno)d - %(message)s")) # 设置格式
    logger.addHandler(ch)

    # this prevents double outputs
    logger.propagate = False
    return logger
    

logger = _setup_logger()

# 修改logger中handler的日志格式
def format_logger(format):
    for handler in logger.handlers:
        handler.setFormatter(logging.Formatter(format))



if __name__ == "__main__":
    logger.info("Info Message")
    logger.debug("Debug Message")
    # 尝试修改更加详细的日志格式
    format_logger("%(name)s -%(filename)s:%(lineno)d - %(message)s")
    logger.info("Info Message")
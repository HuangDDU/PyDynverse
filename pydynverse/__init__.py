from ._settings import settings
from ._logging import logger
from . import wrap
from . import methods
from . import plot
from . import dimred

__all__ = [
    "wrap",
    "methods",
    "plot",
    "dimred",
]

# ASCII艺术字链接：http://patorjk.com/software/taag
logo = """
  _____       _____                                     
 |  __ \     |  __ \                                    
 | |__) |   _| |  | |_   _ _ ____   _____ _ __ ___  ___ 
 |  ___/ | | | |  | | | | | '_ \ \ / / _ \ '__/ __|/ _ \\
 | |   | |_| | |__| | |_| | | | \ V /  __/ |  \__ \  __/
 |_|    \__, |_____/ \__, |_| |_|\_/ \___|_|  |___/\___|
         __/ |        __/ |                             
        |___/        |___/         
"""

print(logo)
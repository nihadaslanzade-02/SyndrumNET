"""
Structured logging for SyndrumNET pipeline.

Provides consistent logging across all modules with automatic file/console output.
"""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


def setup_logger(
    name: str = "syndrumnet",
    log_dir: Optional[Path] = None,
    level: int = logging.INFO,
    console: bool = True,
) -> logging.Logger:
    """
    Setup structured logger with file and console handlers.
    
    Parameters
    ----------
    name : str
        Logger name.
    log_dir : Path, optional
        Directory for log files. If None, logs to 'logs/'.
    level : int
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    console : bool
        Whether to log to console in addition to file.
        
    Returns
    -------
    logging.Logger
        Configured logger instance.
        
    Examples
    --------
    >>> logger = setup_logger('syndrumnet')
    >>> logger.info("Starting pipeline")
    >>> logger.warning("Missing data for 5 genes")
    >>> logger.error("Network construction failed", exc_info=True)
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Create formatters
    detailed_formatter = logging.Formatter(
        fmt='%(asctime)s | %(name)s | %(levelname)s | %(funcName)s:%(lineno)d | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    simple_formatter = logging.Formatter(
        fmt='%(asctime)s | %(levelname)s | %(message)s',
        datefmt='%H:%M:%S'
    )
    
    # File handler
    if log_dir is None:
        log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"syndrumnet_{timestamp}.log"
    
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)  # Always log everything to file
    file_handler.setFormatter(detailed_formatter)
    logger.addHandler(file_handler)
    
    # Console handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(simple_formatter)
        logger.addHandler(console_handler)
    
    # Log setup completion
    logger.info(f"Logger initialized: {name}")
    logger.info(f"Log file: {log_file}")
    
    return logger


class LoggerMixin:
    """
    Mixin to add logger attribute to classes.
    
    Usage
    -----
    class MyClass(LoggerMixin):
        def __init__(self):
            self.logger.info("Initialized")
    """
    
    @property
    def logger(self) -> logging.Logger:
        """Get or create logger for this class."""
        name = f"syndrumnet.{self.__class__.__module__}.{self.__class__.__name__}"
        return logging.getLogger(name)
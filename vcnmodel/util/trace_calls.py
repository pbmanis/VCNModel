import functools
import sys
import time
from functools import wraps

from rich.text import Text
from rich.console import Console
console = Console()
"""
A simple module to allow tracking and timing of calls using decorators.

"""

class TraceCalls(object):
    """ Use as a decorator on functions that should be traced. Several
        functions can be decorated - they will all be indented according
        to their call depth.
        from: https://eli.thegreenplace.net/2012/08/22/easy-tracing-of-nested-function-calls-in-python
    """

    def __init__(
        self, stream=sys.stdout, indent_step=2, show_args=False, show_ret=False,
                show=True):
        self.stream = stream
        self.indent_step = indent_step
        self.show_ret = show_ret
        self.show_args = show_args
        self.show = show

        # This is a class attribute since we want to share the indentation
        # level between different traced functions, in case they call
        # each other.
        TraceCalls.cur_indent = 0

    def __call__(self, fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            indent = " " * TraceCalls.cur_indent
            if self.show_args:
                argstr = ", ".join(
                    [repr(a) for a in args]
                    + ["%s=%s" % (a, repr(b)) for a, b in kwargs.items()]
                )
                if self.show:
                    self.stream.write("-->%s%s(%s)\n" % (indent, fn.__name__, argstr))
            else:
                # self.stream.write('-->%s%s\n' % (indent, fn.__name__))
                text = Text.assemble(
                    (f"* {indent:s}-->{fn.__name__:s}\n", "bold yellow")
                )
                if self.show:
                    console.print(text)

            TraceCalls.cur_indent += self.indent_step
            ret = fn(*args, **kwargs)
            TraceCalls.cur_indent -= self.indent_step

            if self.show_ret and self.show:
                self.stream.write("%s returns--> %s\n" % (indent, ret))
            return ret

        return wrapper


def winprint(func):
    """
    Wrapper decorator for functions that print to the text area
    Clears the print area first,
    and puts a line of '*' when the function returns
    """

    @functools.wraps(func)
    def wrapper_print(self, *args, **kwargs):
        self.textclear()
        value = func(self, *args, **kwargs)
        # end_time = time.perf_counter()      # 2
        # run_time = end_time - start_time    # 3
        self.textappend("*" * 80)
        # print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_print


def winprint_continuous(func):
    """
    Wrapper decorator for functions that print to the text area
    DOES NOT clear the print area first,
    """

    @functools.wraps(func)
    def wrapper_print(self, *args, **kwargs):
        value = func(self, *args, **kwargs)
        # end_time = time.perf_counter()      # 2
        # run_time = end_time - start_time    # 3
        # print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_print


def time_func(func):
    """
    Decorator to ime functions.
    Place inside (after) winprint if using
    Output is to terminal.
    """

    @functools.wraps(func)
    def wrapper_timer(self, *args, **kwargs):
        print(f"Starting : {func.__name__!r}")
        start_time = time.perf_counter()  # 1
        value = func(self, *args, **kwargs)
        end_time = time.perf_counter()  # 2
        run_time = end_time - start_time  # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_timer
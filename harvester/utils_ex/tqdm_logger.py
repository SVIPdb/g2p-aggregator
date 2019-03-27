import tqdm
import contextlib
import sys
import logging


class DummyTqdmFile(object):
    """Dummy file-like that will write to tqdm"""
    file = None

    def __init__(self, f):
        self.file = f

    def write(self, x):
        # Avoid print() second call (useless \n)
        if len(x.rstrip()) > 0:
            tqdm.tqdm.write(x, file=self.file)

    def flush(self):
        return getattr(self.file, "flush", lambda: None)()


@contextlib.contextmanager
def std_out_err_redirect_tqdm():
    """
    Captures anything written to sys.stdout, sys.stderr and redirects it to tqdm.write(). Restores the previous file
    handle assocations when this context manager is exited.

    Borrowed from https://github.com/tqdm/tqdm/blob/master/examples/redirect_print.py, with some extra comments added.
    """
    orig_out_err = sys.stdout, sys.stderr
    try:
        # replaces both streams with instances of DummyTqdmFile applied to each stream
        sys.stdout, sys.stderr = map(DummyTqdmFile, orig_out_err)
        yield orig_out_err[0]
    except Exception as exc:
        # Relay exceptions
        # FIXME: uh, isn't this a vacuous statement?
        raise exc
    finally:
        # Always restore sys.stdout/err if necessary
        sys.stdout, sys.stderr = orig_out_err


class TqdmLoggingHandler(logging.Handler):
    """
    Redirects logging attempts to tqdm.write(), if a tqdm progress bar is currently being displayed.

    Presumably borrowed from https://github.com/tqdm/tqdm/issues/313#issuecomment-347960988
    """
    def __init__(self, level=logging.NOTSET):
        super(TqdmLoggingHandler, self).__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

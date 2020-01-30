# coding=utf-8
import logging
import timeit
from collections import defaultdict

runtime_stats = defaultdict(int)


def show_runtime_stats(bar_width=25):
    print "Runtime stats:"

    if len(runtime_stats) > 0:
        total = sum(runtime_stats.values())
        highest = max(runtime_stats.values())
        max_name_width = max(len(x) for x in runtime_stats.keys())

        for name, elapsed in runtime_stats.items():
            pct = (elapsed / total) * 100.0
            filled_bars = int(elapsed / total * bar_width)
            print u"- %*s: %15.3f: %-*s (%.1f%%)" % (
                max_name_width, name, elapsed, bar_width + 1,
                (u"\u2593" * filled_bars) + (u"\u2591" * (bar_width - filled_bars)),
                pct
            )

    else:
        print "- (no runtime stats collected)"

class OpTimeLogger:
    """
    Emits how long it took to run once it exits
    """
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start_time = timeit.default_timer()
        return self

    def __exit__(self, *exc_args):
        elapsed = max(timeit.default_timer() - self.start_time, 0)
        logging.info("%s took %.3fs to run" % (self.name, elapsed))


class DelayedOpLogger:
    """
    Emits a log statement if the contained code takes longer than 'duration' seconds.

    If 'duration' is none, no logging statement is emitted.
    If 'aggregate_times' is True, also sums aggregate per-normalizer runtime to the 'runtime_stats' global dict.
    """
    def __init__(self, name, duration=1.0, aggregate_times=True):
        self.duration = duration
        self.more = None
        self.name = name
        self.aggregate_times = aggregate_times

    def __enter__(self):
        self.start_time = timeit.default_timer()
        return self

    def logdelayed(self, *more):
        """
        Customizes the log statement with extra information.
        :param more: values to append to the log statement
        """
        self.more = more

    def __exit__(self, *exc_args):
        elapsed = max(timeit.default_timer() - self.start_time, 0)

        if self.aggregate_times:
            runtime_stats[self.name] += elapsed

        if self.duration and elapsed > self.duration:
            if self.more:
                logging.info("%s delayed %.3fs; %s" % (self.name, elapsed, " ".join(str(x) for x in self.more)))
            else:
                logging.info("%s delayed %.3fs" % (self.name, elapsed))


def add_crawl_status(obj, subkey, entry):
    """
    Ensures first that obj has a dict called 'crawl_status', then that it has an array obj['crawl_status'][subkey],
    then concatenates entry to that array.
    :param obj: the object to annotate with a status entry
    :param subkey: the subkey within the status record, e.g. 'location_normalizer' for errors encountered in that module.
    :param entry: the actual log entry to store, e.g. "HGVS parsing failed"
    :return:
    """

    if 'crawl_status' not in obj or obj['crawl_status'] is None:
        obj['crawl_status'] = defaultdict(list)

    obj['crawl_status'][subkey].append(entry)

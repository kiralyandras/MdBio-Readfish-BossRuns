"""_alignment.py

Classes related to alignment
"""
import logging

from ru.utils import setup_logger

logger = setup_logger(__name__, "%(asctime)s %(name)s %(message)s", level=logging.INFO)
_mappy_rs = False
try:
    import mappy_rs as mappy

    _mappy_rs = True
    logger.debug("Using mappy-rs")
except ImportError:
    mappy = None

if mappy is None:
    try:
        import mappy

        logger.debug("Using mappy")
    except ImportError:
        raise ImportError(
            "Cannot find either `mappy-rs` nor `mappy`. One of these is required."
        )


class Mapper:
    def __init__(self, index="", **kwargs):
        logger.debug(f"Opening index: {index!r}")
        self.index = index
        self.aligner = mappy.Aligner(self.index, **kwargs)
        logger.debug(f"Index initialised? {self.initialised}")
        if _mappy_rs:
            self.enable_threading(1)

    @property
    def has_multithreading(self):
        return _mappy_rs

    def enable_threading(self, threads):
        if _mappy_rs:
            self.aligner.enable_threading(threads)
        else:
            raise RuntimeError("Not using mappy-rs, can't use multithreading")

    @property
    def initialised(self):
        return bool(self.aligner)

    def map_reads(self, *args, **kwargs):
        if _mappy_rs:
            yield from self._rust_mappy_wrapper(*args, **kwargs)
        else:
            yield from self._c_mappy_wrapper(*args, **kwargs)

    def _c_mappy_wrapper(self, basecalls, *, key=lambda x: x):
        n = 0
        for meta, data in basecalls:
            yield meta, data, list(self.aligner.map(key(data)))
            n += 1
        logger.debug(f"Got {n:,} alignments")

    def _rust_mappy_wrapper(self, basecalls, *, key=lambda x: x):
        def _gen(_basecalls, _key):
            n_, m = 0, 0
            for meta, data in _basecalls:
                fa = _key(data)
                if not fa:
                    m += 1
                    continue
                yield {"seq": fa, "meta": meta, "basecall": data}
                n_ += 1
            logger.debug(f"Gen sent {n_:,} sequences, skipped {m:,}")

        recv = self.aligner.map_batch(_gen(basecalls, key))
        n = 0
        for mappings, sent_data in recv:
            yield sent_data["meta"], sent_data["basecall"], mappings
            n += 1
        logger.debug(f"Got back {n:,} alignments")

    @property
    def mapper(self):
        return self.aligner

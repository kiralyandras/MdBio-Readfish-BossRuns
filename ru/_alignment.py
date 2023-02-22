"""_alignment.py

Classes related to alignment
"""
import functools

_mappy_rs = False
try:
    import mappy_rs as mappy

    _mappy_rs = True
except ImportError:
    mappy = None

if mappy is None:
    try:
        import mappy
    except ImportError:
        raise ImportError(
            "Cannot find either `mappy-rs` nor `mappy`. One of these is required."
        )


class Mapper:
    def __init__(self, index="", **kwargs):
        self.index = index
        self.aligner = mappy.Aligner(self.index, **kwargs)
        if _mappy_rs:
            self.aligner.enable_threading(1)

    @property
    def initialised(self):
        return bool(self.aligner)

    def map_reads(self, *args, **kwargs):
        if _mappy_rs:
            yield from self._rust_mappy_wrapper(*args, **kwargs)
        else:
            yield from self._c_mappy_wrapper(*args, **kwargs)

    def _c_mappy_wrapper(self, basecalls, *, key=lambda x: x):
        for meta, data in basecalls:
            yield meta, data, list(self.aligner.map(key(data)))

    def _rust_mappy_wrapper(self, basecalls, *, key=lambda x: x):
        def _gen(_basecalls, _key):
            for meta, data in _basecalls:
                yield {"seq": _key(data), "meta": meta, "basecall": data}

        recv = self.aligner.map_batch(_gen(basecalls, key))
        for mappings, sent_data in recv:
            yield sent_data["meta"], sent_data["basecall"], mappings

    @property
    def mapper(self):
        return self.aligner




# f = functools.partial(deep_get, key="datasets.sequence", default="")

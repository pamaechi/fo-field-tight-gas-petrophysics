"""Repair TIF-wrapped LIS tapes whose TIF records carry padding beyond the LIS physical-record length.
Strips the 12-byte TIF headers, trims every physical record to the length declared in its LIS PRH,
and writes a plain LIS byte stream that dlisio reads cleanly.  Usage: tif_repair(src, dst) -> stats dict"""
import struct

def tif_repair(src, dst):
    data = open(src, "rb").read()
    out = bytearray()
    pos, n, padded, marks = 0, 0, 0, 0
    while pos + 12 <= len(data):
        typ, prev, nxt = struct.unpack("<III", data[pos:pos + 12])
        if nxt <= pos or nxt > len(data):          # end of data / corrupt pointer
            break
        rec = data[pos + 12:nxt]
        pos = nxt
        if typ == 1 or len(rec) < 4:               # tape mark or empty
            marks += 1
            continue
        prh_len = struct.unpack(">H", rec[:2])[0]
        if 4 <= prh_len <= len(rec):
            if prh_len < len(rec):
                padded += 1
            rec = rec[:prh_len]
        out += rec
        n += 1
    open(dst, "wb").write(out)
    return dict(records=n, padded=padded, tape_marks=marks, bytes_in=len(data), bytes_out=len(out))

if __name__ == "__main__":
    import sys
    print(tif_repair(sys.argv[1], sys.argv[2]))

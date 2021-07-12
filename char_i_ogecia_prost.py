def osi_mom_bezw_prost(b, h, off=0):
    return b * (h ** 3) /12 + b * h * (off ** 2)

mom_bezw = osi_mom_bezw_prost(0.2, 0.3, 0)

def ugiecie_belki_swob_q(l, q, modyounga, mom_bezw):
    return 5/384 * q * (l ** 4) / (modyounga * mom_bezw)

print(ugiecie_belki_swob_q(5, 0.007, 31, mom_bezw))
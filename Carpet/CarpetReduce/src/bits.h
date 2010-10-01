#ifndef BITS_H
#define BITS_H

/* a mask with exactly bit n set */
#define BMSK(n)   (1U << (n))

/* whether bit n in value v is set */
#define BGET(v,n) (!!((v) & BMSK(n)))

/* set, clear, or invert bit n in value v */
#define BSET(v,n) ((v) |  BMSK(n))
#define BCLR(v,n) ((v) & ~BMSK(n))
#define BINV(v,n) ((v) ^  BMSK(n))

/* set bit n in value v to b */
#define BCPY(v,b,n) ((b) ? BSET(v,n) : BCLR(v,n))

/* count the number of set bits */
#define BCNT(x) (((BX_(x)+(BX_(x)>>4)) & 0x0F0F0F0F) % 0xFF)
#define BX_(x)  ((x)                                            \
                 - (((x)>>1)&0x77777777)                        \
                 - (((x)>>2)&0x33333333)                        \
                 - (((x)>>3)&0x11111111))

#endif  /* #ifndef BITS_H */

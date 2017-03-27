#ifdef NELEMENTS
#undef NELEMENTS
#endif

#define NELEMENTS(x) (sizeof(x)/sizeof(x[0]))

#ifdef _DEBUG
#define MYASSERT_DEFINED 1
#endif

#if MYASSERT_DEFINED
inline void MyAssert(unsigned condition)
{
    if (!condition)
    {
        __debugbreak();
    }
}
#else
#define MyAssert(x)
#endif

// make n steps to the right from f
inline float next_float(float f, unsigned n = 1)
{
    float fNext;
    if (f > 0)
    {
        (unsigned &)fNext = (unsigned &)f + n;
    }
    else
    {
        (unsigned &)fNext = (unsigned &)f - n;
    }
    MyAssert(fNext > f);
    return fNext;
}
// make n steps to the left from f
inline float prev_float(float f, unsigned n = 1)
{
    float fPrev;
    if (f > 0)
    {
        (unsigned &)fPrev = (unsigned &)f - n;
    }
    else
    {
        (unsigned &)fPrev = (unsigned &)f + n;
    }
    MyAssert(fPrev < f);
    return fPrev;
}

// min and max from std:: conflict with those defined elsewhere. that happens in multiple places - so we undefine them here to zap that problem once and for all
// just include MyMisc.h everywhere
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#ifdef SAFE_RELEASE
#undef SAFE_RELEASE
#endif

template <class T>
inline void SAFE_RELEASE(T *(&p))
{
  if (p)
  {
    p->Release();
    p = NULL;
  }
}

template <class T, unsigned N>
inline unsigned NUMELEMENTS(T (&)[N])
{
  return N;
}

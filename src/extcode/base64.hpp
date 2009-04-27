class CBase64Encode
{
private:
  static const char encodeCharacterTable[];
  std::ostream& out;

  unsigned short i;
  char buff1[3];

public:
  CBase64Encode(std::ostream& nout):
    out(nout), i(0)
  {}

  ~CBase64Encode()
  {
    char buff2[4];

    if (i)
      {
	for(unsigned short j=i;j<3;j++) buff1[j] = '\0';
	buff2[0] = (buff1[0] & 0xfc) >> 2;
	buff2[1] = ((buff1[0] & 0x03) << 4) + ((buff1[1] & 0xf0) >> 4);
	buff2[2] = ((buff1[1] & 0x0f) << 2) + ((buff1[2] & 0xc0) >> 6);
	buff2[3] = buff1[2] & 0x3f;
	for (unsigned short j=0;j<(i+1);j++) out << encodeCharacterTable[buff2[j]];
	while(i++<3) out << '=';
      }
  }

  CBase64Encode& operator<<(const char* val)
  {
    unsigned short j(0);
    while((buff1[i++] = val[j++]) != '\0')
      if (i==3)
	{
	  out << encodeCharacterTable[(buff1[0] & 0xfc) >> 2];
	  out << encodeCharacterTable[((buff1[0] & 0x03) << 4) + ((buff1[1] & 0xf0) >> 4)];
	  out << encodeCharacterTable[((buff1[1] & 0x0f) << 2) + ((buff1[2] & 0xc0) >> 6)];
	  out << encodeCharacterTable[buff1[2] & 0x3f];
	  i=0;
	}

    --i;
    return *this;
  }

  template<class T>
  CBase64Encode& operator<<(const T& val)
  {
    for (unsigned short j=0; j < sizeof(T); ++j)
      {
	buff1[i++] = ((const char*)&val)[j];
	if (i==3)
	  {
	    out << encodeCharacterTable[(buff1[0] & 0xfc) >> 2];
	    out << encodeCharacterTable[((buff1[0] & 0x03) << 4) + ((buff1[1] & 0xf0) >> 4)];
	    out << encodeCharacterTable[((buff1[1] & 0x0f) << 2) + ((buff1[2] & 0xc0) >> 6)];
	    out << encodeCharacterTable[buff1[2] & 0x3f];
	    i=0;
	  }
      }

    return *this;
  }

};

class CBase64Decode
{
private:
  static const char decodeCharacterTable[];
  std::istream& in;

public:
  CBase64Decode(std::istream& nin):
    in(nin) {}

  void Decode(std::ostream &out)
  {
    char buff1[4];
    char buff2[4];
    unsigned short i=0, j;
    while(in.readsome(&buff2[i], 1) && buff2[i] != '=')
      {
	if (++i==4)
	  {
	    for (i=0;i!=4;i++)
	      buff2[i] = decodeCharacterTable[buff2[i]];
	    out << (char)((buff2[0] << 2) + ((buff2[1] & 0x30) >> 4));
	    out << (char)(((buff2[1] & 0xf) << 4) + ((buff2[2] & 0x3c) >> 2));
	    out << (char)(((buff2[2] & 0x3) << 6) + buff2[3]);
	    i=0;
	  }
      }
    if (i)
      {
	for (j=i;j<4;j++) buff2[j] = '\0';
	for (j=0;j<4;j++) buff2[j] = decodeCharacterTable[buff2[j]];
	buff1[0] = (buff2[0] << 2) + ((buff2[1] & 0x30) >> 4);
	buff1[1] = ((buff2[1] & 0xf) << 4) + ((buff2[2] & 0x3c) >> 2);
	buff1[2] = ((buff2[2] & 0x3) << 6) + buff2[3];
	for (j=0;j<(i-1); j++) out << (char)buff1[j];
      }
  }
};

#ifndef __RANDGEN_H__
#define __RANDGEN_H__

class RandGen {
	public:
		RandGen(){}
		virtual ~RandGen(){}
		virtual double randomOne() = 0;
	};

#endif

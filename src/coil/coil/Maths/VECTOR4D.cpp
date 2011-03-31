//////////////////////////////////////////////////////////////////////////////////////////
//	VECTOR4D.cpp
//	Function definitions for 4d vector class
//	Downloaded from: www.paulsprojects.net
//	Created:	20th July 2002
//	Modified:	15th August 2002	-	prevent divide by zero in operator VECTOR3D()
//	Modified:	8th November 2002	-	Changed Constructor layout
//									-	Some speed Improvements
//									-	Corrected Lerp
//
//	Copyright (c) 2006, Paul Baker
//	Distributed under the New BSD Licence. (See accompanying file License.txt or copy at
//	http://www.paulsprojects.net/NewBSDLicense.txt)
//////////////////////////////////////////////////////////////////////////////////////////	

#include "Maths.h"

VECTOR4D operator*(float scaleFactor, const VECTOR4D & rhs)
{
	return rhs*scaleFactor;
}

bool VECTOR4D::operator==(const VECTOR4D & rhs) const
{
	if(x==rhs.x && y==rhs.y && z==rhs.z && w==rhs.w)
		return true;

	return false;
}

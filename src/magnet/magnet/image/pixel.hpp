/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Severin Strobl <-->

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <inttypes.h>

namespace magnet {
  namespace image {
    enum ColorType { RGB = PNG_COLOR_TYPE_RGB, RGBA = PNG_COLOR_TYPE_RGBA };

    template<ColorType = RGB>
    class Pixel;

    template<>
    class Pixel<RGB> {
    public:
      Pixel() : r_(0), g_(0), b_(0) {
      }

      Pixel(uint8_t r, uint8_t g, uint8_t b) : r_(r), g_(g), b_(b) {
      }

      uint8_t red() const {
	return r_;
      }

      uint8_t& red() {
	return r_;
      }

      void red(const uint8_t r) {
	r_ = r;
      }

      uint8_t green() const {
	return g_;
      }

      uint8_t& green() {
	return g_;
      }

      void green(const uint8_t g) {
	g_ = g;
      }

      uint8_t blue() const {
	return b_;
      }

      uint8_t& blue() {
	return b_;
      }

      void blue(const uint8_t b) {
	b_ = b;
      }

      void convert(const png_bytep pngPixel) {
	r_ = pngPixel[0];
	g_ = pngPixel[1];
	b_ = pngPixel[2];
      }

    private:
      uint8_t r_, g_, b_;
    };

    template<>
    class Pixel<RGBA> {
    public:
      Pixel() : r_(0), g_(0), b_(0), a_(0) {
      }

      Pixel(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : r_(r), g_(g),
							  b_(b), a_(a) {
      }

      uint8_t red() const {
	return r_;
      }

      uint8_t& red() {
	return r_;
      }

      void red(const uint8_t r) {
	r_ = r;
      }

      uint8_t green() const {
	return g_;
      }

      uint8_t& green() {
	return g_;
      }

      void green(const uint8_t g) {
	g_ = g;
      }

      uint8_t blue() const {
	return b_;
      }

      uint8_t& blue() {
	return b_;
      }

      void blue(const uint8_t b) {
	b_ = b;
      }

      uint8_t alpha() const {
	return a_;
      }

      uint8_t& alpha() {
	return a_;
      }

      void alpha(const uint8_t a) {
	a_ = a;
      }

      void convert(const png_bytep pngPixel) {
	r_ = pngPixel[0];
	g_ = pngPixel[1];
	b_ = pngPixel[2];
	a_ = pngPixel[3];
      }

    private:
      uint8_t r_, g_, b_, a_;
    };
  }
}

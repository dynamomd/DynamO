/* Compile this with the following code
   g++ ffmpeg-example.c  $(pkg-config  --cflags libavcodec) $(pkg-config --libs libavcodec) -I ~/dynamo/src/magnet/ -O3
*/

#include <magnet/image/videoEncoder.hpp>

int main(int argc, char **argv)
{
  const char *filename = "/tmp/test.mpg";
  size_t width = 1023;
  size_t height = 1023;
  size_t size = width * height;

  VideoEncoder encoder(filename, width, height);

  /* Allocate the RGB image in advance */
  std::vector<uint8_t> rgb_buf(size * 3);
  
  /* encode 3 seconds of video */
  for(int i=0; i<75; ++i) 
    {
      for(int y=0; y < height; y++)
	for(int x=0; x< width; x++)
	  {
	    rgb_buf[3*(x + width * y) + 0] = 10*i; //R
	    rgb_buf[3*(x + width * y) + 1] = 51*i; //G
	    rgb_buf[3*(x + width * y) + 2] = std::min(x,y); //B
	  }
      encoder.addFrame(rgb_buf);
    }
  
  encoder.close();

  return 0;
}

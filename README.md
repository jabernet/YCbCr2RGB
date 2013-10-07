YCbCr to RGB
=============

In conversion.cpp there is a simple implementation that converts seperated inputs of Y Cb and Cr channels from the Theora codec into RGB output.
The conversion will operate on 16 pixels in parallel by loading them as epi8 into an __m128 value. They are then converted to 32bit floats, the operations specified in http://www.theora.org/doc/Theora.pdf applied and assembled as an RGB stream into the output video.
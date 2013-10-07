typedef struct{
  /**The width of this plane.*/
  int            width;
  /**The height of this plane.*/
  int            height;
  /**The offset in bytes between successive rows.*/
  int            stride;
  /**A pointer to the beginning of the first row.*/
  unsigned char *data;
}th_img_plane;

typedef unsigned __int32 ogg_uint32_t;

typedef enum{
  TH_CS_UNSPECIFIED,
  TH_CS_ITU_REC_470M,
  TH_CS_ITU_REC_470BG,
  TH_CS_NSPACES
}th_colorspace;

typedef enum{
  TH_PF_420,
  TH_PF_RSVD,
  TH_PF_422,
  TH_PF_444,
  TH_PF_NFORMATS
}th_pixel_fmt;

typedef struct{
  unsigned char version_major;
  unsigned char version_minor;
  unsigned char version_subminor;
  ogg_uint32_t  frame_width;
  ogg_uint32_t  frame_height;
  ogg_uint32_t  pic_width;
  ogg_uint32_t  pic_height;
  ogg_uint32_t  pic_x;
  ogg_uint32_t  pic_y;
  ogg_uint32_t  fps_numerator;
  ogg_uint32_t  fps_denominator;
  ogg_uint32_t  aspect_numerator;
  ogg_uint32_t  aspect_denominator;
  th_colorspace colorspace;
  th_pixel_fmt  pixel_fmt;
  int           target_bitrate;
  int           quality;
  int           keyframe_granule_shift;
}th_info;

typedef th_img_plane th_ycbcr_buffer[3];
#ffmpeg -framerate 25 -i  GD424_GD424_forvid_p%02d00_compositionmk2.png  -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
ffmpeg -framerate 25 -i  Example_%03d_composition_rel_He.png  -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" out.mp4

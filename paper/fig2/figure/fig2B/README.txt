## Making each frame
foreach i (`seq 0 39`)
mkdir run/run$i
python render.py ../../simulation/wt/out/run${i}_control_dX0_0_tstart_0.csv run/run$i
end

## Making movie
foreach i (`seq 0 39`)
ffmpeg -framerate 30 -i ./run/run${i}/frame_%04d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p ./run/run${i}/run${i}.mp4

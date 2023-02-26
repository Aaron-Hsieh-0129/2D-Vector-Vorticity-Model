cd src && python PlotAll.py

cd ../graphs/qc+qr+th+u+w
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p qc+qr+th+u+w.mov -y
ffmpeg -i qc+qr+th+u+w.mov -pix_fmt rgb24 qc+qr+th+u+w.gif -y

cd ../zeta
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p zeta.mov -y
ffmpeg -i zeta.mov -pix_fmt rgb24 zeta.gif -y

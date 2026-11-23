"""Plot the signed Benchmark B misfit changes computed by analyze.py."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

out=Path(__file__).resolve().parent
r=json.loads((out/"residual_summary.json").read_text())
b=r["channels"]["fused"]["frequency_bands"]
fig, (ax, bx)=plt.subplots(1,2,figsize=(12,4.7),layout="constrained")
x=np.arange(len(b))
for offset,key,label,color in [(-.18,"delta_ncc","NCC","#166d83"),(.18,"delta_w2","Wasserstein","#bd641c")]:
    ax.bar(x+offset,[i[key] for i in b],width=.35,label=label,color=color)
ax.axhline(0,color="#333333",lw=.8)
ax.set_xticks(x,[i["band"] for i in b],rotation=22,ha="right")
ax.set(ylabel="Kontribusi perubahan misfit terhadap prior",title="Perubahan misfit per pita frekuensi")
ax.legend(frameon=False)
comp=r["components"];x=np.arange(len(comp))
for off,key,label,color in [(-.25,"center_J","Pusat komponen","#4c638a"),
                           (0,"conditional_mean_spectrum_J","Prediksi per komponen","#258375"),
                           (.25,"leave_one_component_out_J","Campuran tanpa komponen ini","#ae6c32")]:
    bx.bar(x+off,[i[key] for i in comp],width=.23,label=label,color=color)
for key,label,style in [("prior","Model tengah prior","--"),("prediction","Prediksi campuran penuh",":")]:
    bx.axhline(r["channels"]["fused"][key],ls=style,lw=1.4,color="#143047",label=label)
bx.set_xticks(x,[str(i["component"]) for i in comp])
bx.set_ylim(0, .77)
bx.set(xlabel="Komponen campuran",ylabel="Misfit spektrum gabungan",title="Kecocokan komponen dan campuran")
bx.legend(frameon=False,fontsize=8)
fig.suptitle("Benchmark B | Diagnostik hasil v6",fontsize=15)
fig.savefig(out/"diagnostik_residual_B.png",dpi=170)

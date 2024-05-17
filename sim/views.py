from django.shortcuts import render

from django.http import HttpResponse, HttpRequest, HttpResponseBadRequest, HttpResponseServerError
from yee3d.Grid import Grid
from yee3d.BoundaryCond import FirstOrderABC
from django.views.decorators.csrf import csrf_exempt
import json
import io
import base64

from yee3d.GridVisualizer import GridVisualizer
from PIL import Image

# Create your views here.
@csrf_exempt
def index(request: HttpRequest):
    if request.method == "GET":
        return HttpResponse("Hi")
    if request.method == "POST":
        config = json.loads(request.body)
        if config == None:
            return HttpResponseBadRequest("Config not found.")
        
        print(config)
        G = Grid(int(config['X']), int(config['Y']), int(config['Z']), int(config['maxTime']))
        abc = FirstOrderABC(G)

        V = GridVisualizer(G, 'test')
        path = V.snapshotGIF(G.maxTime, period=1, dim=0, slice_index=(G.dims[0] - 1)//2, field='Hx')
        # print(path)

        try:
            with open(path, 'rb') as f:
                gif_file = f.read()
                detect_base64 = 'data:image/gif;base64,{}'.format(base64.b64encode(gif_file).decode())
                V.cleanup()
                return HttpResponse(detect_base64)
        except IOError:
            return HttpResponseServerError("An error occurred generating the simulation gif.")



from django.shortcuts import render

from django.http import HttpResponse, HttpRequest, HttpResponseBadRequest, HttpResponseServerError
from yee3d.Grid import Grid
from yee3d.BoundaryCond import FirstOrderABC
from django.views.decorators.csrf import csrf_exempt
import json

from yee3d.GridVisualizer import GridVisualizer
from PIL import Image

# Create your views here.
@csrf_exempt
def index(request: HttpRequest):
    if request.method == "GET":
        return HttpResponse("Hi")
    if request.method == "POST":
        body = json.loads(request.body)
        config = body.get('gridConfig', None)
        if config == None:
            return HttpResponseBadRequest("Config not found.")
        
        # print(config)
        
        G = Grid(int(config['sx']), int(config['sy']), int(config['sz']), int(config['maxTime']))
        abc = FirstOrderABC(G)

        V = GridVisualizer(G, 'test')
        path = V.snapshotGIF(G.maxTime, period=1, dim=0, slice_index=(G.dims[0] - 1)//2, field='Hx')
        # print(path)

        try:
            with open(path, "rb") as f:
                r = f.read()
                V.cleanup()
                return HttpResponse(r, content_type="image/gif")
        except IOError:
            return HttpResponseServerError("An error occurred generating the simulation gif.")



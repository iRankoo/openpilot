#!/usr/bin/python
#Created by iRankoo

import threading
import requests
from http.server import BaseHTTPRequestHandler, HTTPServer
from os import curdir, sep
import cgi
import json
from openpilot.common.params import Params

PORT_NUMBER = 8888

#This class will handles any incoming request from
#the browser
class myHandler(BaseHTTPRequestHandler):
    #Handler for the GET requests
    def do_GET(self):
        print(self.path)
        if self.path=="/":
            self.path="./nav.html"
        elif self.path=="/nav":
            self.path="./nav.html"

        try:
            #Check the file extension required and
            #set the right mime type

            sendReply = False
            if self.path.endswith(".html"):
                mimetype='text/html'
                sendReply = True
            if self.path.endswith(".jpg"):
                mimetype='image/jpg'
                sendReply = True
            if self.path.endswith(".gif"):
                mimetype='image/gif'
                sendReply = True
            if self.path.endswith(".js"):
                mimetype='application/javascript'
                sendReply = True
            if self.path.endswith(".css"):
                mimetype='text/css'
                sendReply = True

            if sendReply == True:
                #Open the static file requested and send it
                f = open(curdir + sep + self.path)
                print(curdir)
                print(sep)
                print(self.path)

                self.send_response(200)
                self.send_header('Content-type',mimetype)
                self.end_headers()
                buf = f.read()
                byte_data = buf.encode('utf-8')
                self.wfile.write(byte_data)
                f.close()
            return

        except IOError:
            self.send_error(404,'File Not Found: %s' % self.path)

    #Handler for the POST requests
    def do_POST(self):
        if self.path=="/nav":
            form = cgi.FieldStorage(
                fp=self.rfile,
                headers=self.headers,
                environ={'REQUEST_METHOD':'POST',
                         'CONTENT_TYPE':self.headers['Content-Type'],
            })


            try:
                print ("Map url is: %s" % form["nav_url"].value)
            except:
                self.send_response(200)
                self.end_headers()
                buf = "Pls input the map url!"
                byte_data = buf.encode('utf-8')
                self.wfile.write(byte_data)
                return

            map_rul = form["nav_url"].value
            lat, lng = getLatNLng(map_rul)
            if lat == None:
                url = getRedirectUrl(map_rul)
                lat, lng = getLatNLng_Amap(url)

            if lat==None or lng==None:
                self.send_response(200)
                self.end_headers()
                buf = "Pls input the map url!"
                byte_data = buf.encode('utf-8')
                self.wfile.write(byte_data)
                return
            else:
                writeParam(lat, lng)

            self.send_response(200)
            self.end_headers()
            buf = "Thank you and enjoy the c3!"
            byte_data = buf.encode('utf-8')
            self.wfile.write(byte_data)
            return


def getLatNLng(map_rul):
    lat = None
    lng = None
    match = "place?ll="
    iFrom = map_rul.find(match)
    if iFrom > 0:
        iEnd = map_rul.find("&q=", iFrom+len(match))
        ll = map_rul[iFrom+len(match):iEnd]

        c = ll.find(",")
        if c > 0:
            lat = ll[:c]
            lng = ll[c+1:]
            print(lat+" "+lng)
    return lat, lng

def getLatNLng_Amap(map_rul):
    lat = None
    lng = None
    match = "amap.com/?q="
    iFrom = map_rul.find(match)
    if iFrom > 0:
        iEnd = map_rul.find(",%", iFrom+len(match))
        ll = map_rul[iFrom+len(match):iEnd]

        c = ll.find(",")
        if c > 0:
            lat = ll[:c]
            lng = ll[c+1:]
            print(lat+" "+lng)
    return lat, lng

def getRedirectUrl(url):
    r = requests.get(url,headers={"Content-Type":"application/json"})
    reditList = r.history
    print(f'First : {reditList}')
    print(f'Redirect Header : {reditList[0].headers}')
    fUrl = reditList[len(reditList)-1].headers["location"]
    print(f'Final url : {fUrl}')
    return fUrl

def writeParam(lat,lng):
    data = {
    "save_type": "recent",
    "label": "work",
    "place_name": "",
    "place_details": "",
    "time": 0,
    "latitude": float(lat),
    "longitude": float(lng)
    }

    jsn = json.dumps(data)
    print(jsn)

    params = Params()
    params.put("NavDestination", jsn)

def threadHttp():
    try:
        #Create a web server and define the handler to manage the
        #incoming request
        server = HTTPServer(('', PORT_NUMBER), myHandler)
        print ('Started httpserver on port ' , PORT_NUMBER)

        #Wait forever for incoming htto requests
        server.serve_forever()

    except KeyboardInterrupt:
        print ('^C received, shutting down the web server')
        server.socket.close()

def main():
    t = threading.Thread(target=threadHttp)
    t.start()
    t.join()

if __name__ == "__main__":
  main()

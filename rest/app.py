#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018  <@rubyriver.local>
#
# Distributed under terms of the MIT license.

from flask import Flask
from flask_restful import Resource, Api
import json

app = Flask(__name__)
api = Api(app)

class Model(Resource):
    def get(self):
        read_data = None
        with open('streamflows-example.json', 'r') as f:
            read_data = f.read()
        if read_data:
            return json.dumps(json.loads(read_data))
        # else
        # return 404

api.add_resource(Model, '/model')

if __name__ == '__main__':
    app.run(debug=True)


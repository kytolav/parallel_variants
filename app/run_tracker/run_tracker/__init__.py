# -*- coding: utf-8 -*-

from flask import Flask, url_for, render_template, request, session, escape, redirect

app = Flask(__name__)

from . import views

#
# KEY
#
app.secret_key = '\x17\xa1p\x926\xb5\xad\xbf\x88\x9c\xf0\x7f\xc8\x81e\\5\x9d\xf2M\xc7y\xfb'


if __name__ == "__main__":
    app.run()

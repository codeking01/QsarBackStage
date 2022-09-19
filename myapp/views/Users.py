# auth: code_king
# time: 2022/9/19 14:58
# file: Users.py
from django.shortcuts import render


# 后台的首页
def index(request):
    return render(request, 'index.html')

# auth: code_king
# time: 2022/8/29 19:35
# file: startPages.py
from io import BytesIO

from django.http import HttpResponse
from django.shortcuts import render, redirect

from myapp.models import TestResult
from myapp.utils.SaturatedVaporPressure.predict_Ps import cal_Ps
from myapp.utils.form import PredictModelForm


def predict_test(request):
    if request.method == 'GET':
        # 获取所有城市数据
        queryset = TestResult.objects.all()
        form = PredictModelForm()
        context = {
            'queryset': queryset,
            'form': form
        }
        return render(request, 'predictModel.html', context)
    else:
        template = float(request.POST.get('template'))
        boiling_point = float(request.POST.get('boiling-point'))
        GjfFiles = request.FILES.get('file')
        try:
            f = GjfFiles.readlines()
            f = [item.decode() for item in f]
            GjfFiles.close()
            # 调用事先写好的api，得到预测结果
            P_cul = str(cal_Ps(f, boiling_point, template)).replace('[', '').replace(']', '')
        except Exception as e:
            print(e)
            P_cul= 'Gjf可能存在问题，请联系开发人员！'
        data = {
            'template': template,
            'boiling_point': boiling_point,
            'GjfFiles': GjfFiles.name,
            'results': P_cul
        }
        return render(request, 'predictModel.html', {"data": data})
        # form = PredictModelForm(data=request.POST, files=request.FILES)
        # 拿到数据直接开始渲染
        # if form.is_valid():
        # 对于文件：自动保存；
        # 字段 + 上传路径写入到数据库
        # form.save()
        # 重定向到当前页面，浆液门面数据展示
        # return redirect('/index/')
    # return render(request, 'predictModel.html', {"form": form})
    return render(request, 'predictModel.html')

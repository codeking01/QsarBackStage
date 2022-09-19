# auth: code_king
# time: 2022/8/29 19:35
# file: AllPredModels.py

from django.shortcuts import render

from myapp.utils.Cyclodextrin.Alpha.predict_lgKs import cal_alpha_ks
from myapp.utils.Cyclodextrin.Beta.predict_lgKs import cal_beta_ks
from myapp.utils.SaturatedVaporPressure.predict_Ps import cal_Ps


def predict_test(request):
    if request.method == 'GET':
        # 获取所有城市数据
        # queryset = TestResult.objects.all()
        # form = PredictModelForm()
        # context = {
        #     'queryset': queryset,
        #     'form': form
        # }
        return render(request, 'predictModel.html')
    else:
        template = float(request.POST.get('template'))
        boiling_point = float(request.POST.get('boiling-point'))
        try:
            GjfFiles = request.FILES.get('file')
            f = GjfFiles.readlines()
            f = [item.decode() for item in f]
            GjfFiles.close()
            # 调用事先写好的api，得到预测结果
            P_cul = str(cal_Ps(f, boiling_point, template)).replace('[', '').replace(']', '')
            # P_cul=format(float(P_cul), '.3f')
        except Exception as e:
            print(e)
            template = 'error'
            boiling_point = 'error'
            GjfFiles.name = 'unknown'
            P_cul = 'Gjf可能存在问题，请联系开发人员！'
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


# 环糊精模型
def cyclodextrin_predict(request):
    if request.method == 'GET':
        return render(request, 'cyclodextrin.html')
    else:
        try:
            GjfFiles = request.FILES.get('file')
            f = GjfFiles.readlines()
            f = [item.decode() for item in f]
            GjfFiles.close()
            # 调用事先写好的api，得到预测结果
            # P_cul = str(cal_alpha_ks(f))
            # 计算 alpha常数
            alpha_ks = str(cal_alpha_ks(f)).replace('[', '').replace(']', '')
            # 计算beta常数
            beta_ks = str(cal_beta_ks(f)).replace('[', '').replace(']', '')
            # P_cul=format(float(P_cul), '.3f')
        except Exception as e:
            print(e)
            GjfFiles.name = 'unknown'
            alpha_ks = '存在位置问题，请联系开发人员！'
            beta_ks = '存在位置问题，请联系开发人员！'
        data = {
            'GjfFiles': GjfFiles.name,
            'alpha_ks': alpha_ks,
            'beta_ks': beta_ks
        }
        return render(request, 'cyclodextrin.html', {"data": data})

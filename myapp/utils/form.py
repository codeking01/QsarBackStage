# auth: code_king
# time: 2022/8/30 9:26
# file: form.py
from django import forms
from django.forms import ModelForm

from myapp.models import TestResult


#  预测结果的Form
class PredictModelForm(ModelForm):
    class Meta:
        model = TestResult
        fields = '__all__'
        widgets = {
            'template': forms.TextInput(attrs={'class': 'form-control', 'placeholder': '请输入预测温度,单位为k'}),
            'files': forms.FileInput(attrs={'class': 'form-control', 'type': 'file', 'data-overwrite-initial': "false",
                                            'data-min-file-count': "1"}),
        }

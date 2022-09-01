from django.db import models


# Create your models here.
class TestResult(models.Model):
    """城市"""
    files = models.FileField(verbose_name='gjf文件', upload_to='GjfFiles')
    template = models.FloatField(verbose_name='预测温度', max_length=32)
    results = models.CharField(verbose_name='预测结果', max_length=32)

# Generated by Django 4.0.6 on 2022-08-30 01:48

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='PredicateResult',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('template', models.FloatField(max_length=32, verbose_name='预测温度')),
                ('files', models.FileField(upload_to='GjfFiles', verbose_name='gjf文件')),
                ('results', models.FloatField(max_length=32, verbose_name='预测结果')),
            ],
        ),
    ]

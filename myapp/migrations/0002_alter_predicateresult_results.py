# Generated by Django 4.0.6 on 2022-08-30 06:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('myapp', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='predicateresult',
            name='results',
            field=models.CharField(max_length=32, verbose_name='预测结果'),
        ),
    ]

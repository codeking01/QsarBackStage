"""QsarBackstage URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path

from myapp.views import AllPredModels, deal_pages, Users

urlpatterns = [
    # path('admin/', admin.site.urls),
    path('', Users.index),
    path('saturated_vapor_pressure/', AllPredModels.predict_test),
    # 牛牛的环糊精模型
    path('Cyclodextrin_predict/', AllPredModels.cyclodextrin_predict)
]

handler404 = deal_pages.page_not_found
handler500 = deal_pages.server_error

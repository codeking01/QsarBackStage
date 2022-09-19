# 建立 python3.10 环境
FROM python:3.9

# 镜像作者
MAINTAINER codeking

# 设置 python 环境变量
ENV PYTHONUNBUFFERED 1

# 设置pip源为国内源
COPY pip.conf /root/.pip/pip.conf

# 在容器内创建mysite文件夹
RUN mkdir -p /var/www/projects/qsarBackstage

# 设置容器内工作目录
WORKDIR /var/www/projects/qsarBackstage

# 将当前目录文件加入到容器工作目录中（. 表示当前宿主机目录）
ADD . /var/www/projects/qsarBackstage

# pip安装依赖
RUN pip install -r requirements.txt

from django.shortcuts import render


# Create your views here.
def page_not_found(request, exception=None):
    """
    404 页面
    :param request:
    :return:
    """
    return render(request, "handler_404.html", status=404)


def server_error(request, exception=None):
    """
    500 页面
    :param request:
    :return:
    """
    return render(request, "handler_500.html", status=500)

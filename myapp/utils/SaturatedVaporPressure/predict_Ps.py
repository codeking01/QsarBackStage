import warnings
from numpy import mat
from numpy.linalg import norm
import numpy as np
import math
import copy


def m_s(n_atom, M_S_A, M_adj):
    M_S = M_S_A.copy()
    for m in range(1, n_atom - 1):
        w_ms = np.where(M_S == m)
        if len(w_ms[0]) == 0:  # 节省时间
            break
        for i in range(0, len(w_ms[0])):
            w_msi = M_adj[w_ms[1][i]]
            for j in w_msi:
                if M_S[w_ms[0][i], j - 1] == 0 and w_ms[0][
                    i] != j - 1:  # M_S[w_ms[0][i],j-1]==0去掉描述过的           w_ms[0][i]!=j-1    去掉自身，保证对角线元素为0
                    M_S[w_ms[0][i], j - 1] = m + 1
                    M_S[j - 1, w_ms[0][i]] = m + 1
    return M_S


def msi_gjf(orimsi_all):
    M_atom = []
    n_atom = 0
    for i_L in orimsi_all:
        if i_L[5:12] == '       ':
            n_atom = n_atom + 1
            i_L = i_L.replace('\n', '')
            i_L = i_L.split()
            i_at = i_L[0]
            M_atom.append(i_at)
    M_S_A = mat(np.zeros((n_atom, n_atom)))  # 构建一个相邻矩阵
    M_bon_ = mat(np.zeros((n_atom, n_atom)))  # 构建一个键的矩阵
    for i_L in orimsi_all:
        i_L = i_L.split()
        if len(i_L) > 2 and i_L[0].isdigit():
            i_L_n = list(map(float, i_L))  # 数字全部转成  float
            i_L_n_i = np.hstack((i_L_n[0], i_L_n[1::2])) - 1  # 取原子序号，以0为第一个原子
            i_L_n_i = i_L_n_i.astype(np.int64)  # 转类型，float 转 int
            M_S_A[i_L_n_i[0], list(i_L_n_i[1::])] = 1  # 相连原子标1
            M_bon_[i_L_n_i[0], list(i_L_n_i[1::])] = i_L_n[2::2]  # 相连的方式
    M_S_A = M_S_A + M_S_A.T  # 构建实对称矩阵,大于1的改为1
    # M_S_L=M_S_A
    # for i in range(np.shape(M_S_L)[0]):
    #     M_S_L[i,i]=np.sum(M_S_L[i,:])*(-1)
    M_bon_ = M_bon_ + M_bon_.T  # 需要修改
    M_adj_H = []
    M_adj_nH = []
    M_adj = []
    M_bon = []
    for j in range(0, n_atom):
        M_adj_ = np.where(M_S_A[j, :] > 0)
        # (array([0, 0, 0], dtype=int64), array([ 1,  3, 30], dtype=int64))
        M_adj.append(list(1 + M_adj_[1]))
        M_bon.append(M_bon_[j, M_adj_[1]].tolist()[0])
        M_adj__ = []  # H
        M_adj_nH_ = []  # 非 H
        for i in M_adj_[1]:
            if M_atom[i] == 'H':
                M_adj__.append(i + 1)
            else:
                M_adj_nH_.append(i + 1)
        M_adj_nH.append(M_adj_nH_)
        M_adj_H.append(M_adj__)
    M_S = m_s(n_atom, M_S_A, M_adj)
    M_aro = mat(np.zeros((n_atom, 1)))
    for k in range(n_atom):
        if len(np.where(mat(M_bon[k]) == 1.5)[1]) > 0:
            M_aro[k, 0] = 1
    M_cyc = mat(np.zeros((n_atom, 1)))
    for i in range(0, n_atom):
        if len(M_adj_nH[i]) > 1:  # 为什么不写大于 1
            for j in range(0, n_atom):
                wl = np.where(M_S[j, (mat(M_adj_nH[i]) - 1).tolist()[0]] < M_S[j, i] + 1)[1]
                if len(wl) > 1:
                    M_cyc[i, 0] = 1
                    break
    MSI = {'n_atom': n_atom, 'M_S_A': M_S_A, 'M_S': M_S, 'M_bon': M_bon, 'M_bon_': M_bon_, 'M_atom': M_atom,
           'M_aro': M_aro, 'M_adj': M_adj, 'M_adj_H': M_adj_H, 'M_adj_nH': M_adj_nH, 'M_cyc': M_cyc}
    return MSI


def m_property(name_atom):  # 性质矩阵
    ele_H = 2.2;
    ele_B = 2.04;
    ele_C = 2.55;
    ele_N = 3.04;
    ele_O = 3.44;
    ele_F = 3.98;
    ele_Si = 1.9;
    ele_P = 2.19;
    ele_S = 2.58;
    ele_Cl = 3.16;
    ele_Br = 2.96;
    ele_I = 2.66;
    ele_As = 2.18;
    wei_H = 1.00794;
    wei_B = 10.811;
    wei_C = 12.0107;
    wei_N = 14.0067;
    wei_O = 15.9994;
    wei_F = 18.9984032;
    wei_Si = 28.0855;
    wei_P = 30.973761;
    wei_S = 32.065;
    wei_Cl = 35.453;
    wei_Br = 79.904;
    wei_I = 126.90447;
    wei_As = 74.92160;
    rud_H = 0.79;
    rud_B = 1.17;
    rud_C = 0.91;
    rud_N = 0.75;
    rud_O = 0.65;
    rud_F = 0.57;
    rud_Si = 1.46;
    rud_P = 1.23;
    rud_S = 1.09;
    rud_Cl = 0.97;
    rud_Br = 1.12;
    rud_I = 1.32;
    rud_As = 1.33;
    oue_H = 1;
    oue_B = 3;
    oue_C = 4;
    oue_N = 5;
    oue_O = 6;
    oue_F = 7;
    oue_Si = 4;
    oue_P = 5;
    oue_S = 6;
    oue_Cl = 7;
    oue_Br = 7;
    oue_I = 7;
    oue_As = 5;
    noe_H = 1;
    noe_B = 2;
    noe_C = 2;
    noe_N = 2;
    noe_O = 2;
    noe_F = 2;
    noe_Si = 3;
    noe_P = 3;
    noe_S = 3;
    noe_Cl = 3;
    noe_Br = 4;
    noe_I = 5;
    noe_As = 4;
    ion_H = 1.35984;
    ion_B = 0.82980;
    ion_C = 1.12603;
    ion_N = 1.45341;
    ion_O = 1.36181;
    ion_F = 1.74228;
    ion_Si = 0.81517;
    ion_P = 1.04867;
    ion_S = 1.03600;
    ion_Cl = 1.29676;
    ion_Br = 1.18138;
    ion_I = 1.04513;
    ion_As = 0.97886;
    # ele：电负性     wei:原子质量   rud:原子半径 单位埃   oue:最外层电子数     noe:电子层数     ion:第一电离能，单位：eV/10，类似标准化
    M_ele = [];
    M_wei = [];
    M_rud = [];
    M_oue = [];
    M_noe = [];
    M_ion = [];
    M_aro = []
    for i_L_3 in name_atom:
        exec(
            'M_ele.append(ele_' + i_L_3 + '); M_wei.append(wei_' + i_L_3 + '); M_rud.append(rud_' + i_L_3 + '); M_oue.append(oue_' + i_L_3 + '); M_noe.append(noe_' + i_L_3 + '); M_ion.append(ion_' + i_L_3 + ')')
    MSI = {'M_ele': mat(M_ele).T, 'M_wei': mat(M_wei).T, 'M_rud': mat(M_rud).T, 'M_oue': mat(M_oue).T,
           'M_noe': mat(M_noe).T, 'M_ion': mat(M_ion).T}
    return MSI


def h_delete(in_H, MSI):
    L_H = len(in_H)
    for i, j in MSI.items():
        if j.shape[0] > L_H:
            j = np.delete(j, in_H, 0)
        if j.shape[1] > L_H:
            j = np.delete(j, in_H, 1)
        MSI[i] = j
    return MSI


def hydrogen2heavyatom(M_adj_H, M_P_all, n_atom):
    M_P_all_ = copy.deepcopy(M_P_all)
    for i in M_P_all_:
        for j in range(0, n_atom):
            if not M_adj_H[j]:
                M_P_all_[i][j] = M_P_all_[i][j] + sum(M_P_all_[i][M_adj_H[j]])
    return M_P_all_


def mat_inv(MSI):
    for i, j in MSI.items():
        MSI[i] = 1. / MSI[i]
        MSI[i][MSI[i] == float('inf')] = 0
    return MSI


def norm_d(M, k):
    if k == 0:
        norm_ = math.sqrt(np.min(np.sum(np.power(M, 2), axis=0)))
    elif k == 1:
        norm_ = math.sqrt(np.power(M, 2).sum(axis=0).sum(axis=1))
    elif k == 2:
        norm_ = np.mean(np.sum(abs(M), axis=0))
    elif k == 3:
        norm_ = np.mean(abs(M).sum(axis=0).sum(axis=1))
    elif k == 4:
        norm_ = np.max(np.sum(abs(M), axis=0))
    elif k == 5:
        norm_ = math.sqrt(np.max(np.sum(np.power(M, 2), axis=0)))
    return norm_


def cal_ni(MDS, MP, NI):
    for i in MDS:
        for j in MP:
            MM_1 = np.multiply(i, (j.dot(j.T)))
            MM_2 = np.multiply(i, abs(j - j.T))
            MM_3 = i.dot(j.dot(j.T))
            MM_4 = i.dot(abs(j - j.T))
            NI = np.c_[NI, mat(
                [norm_d(MM_1, 0), norm_d(MM_1, 1), norm_d(MM_1, 2), norm_d(MM_1, 3), norm_d(MM_1, 4), norm_d(MM_1, 5),
                 norm(MM_1, ord=2), norm(MM_1, ord=1),
                 norm_d(MM_2, 0), norm_d(MM_2, 1), norm_d(MM_2, 2), norm_d(MM_2, 3), norm_d(MM_2, 4), norm_d(MM_2, 5),
                 norm(MM_2, ord=2), norm(MM_2, ord=1),
                 norm_d(MM_3, 0), norm_d(MM_3, 1), norm_d(MM_3, 2), norm_d(MM_3, 3), norm_d(MM_3, 4), norm_d(MM_3, 5),
                 norm(MM_3, ord=2), norm(MM_3, ord=1),
                 norm_d(MM_4, 0), norm_d(MM_4, 1), norm_d(MM_4, 2), norm_d(MM_4, 3), norm_d(MM_4, 4), norm_d(MM_4, 5),
                 norm(MM_4, ord=2), norm(MM_4, ord=1)])]
    return NI


def name_ni(name_MDS, name_MP, name_NI):
    for i in name_MDS:
        for j in name_MP:
            MM_1 = i + '.*(' + j + '*' + j + '.T)'
            MM_2 = i + '.*|' + j + '-' + j + '.T|'
            MM_3 = i + '*(' + j + '*' + j + '.T)'
            MM_4 = i + '*|' + j + '-' + j + '.T|'
            name_NI = np.c_[name_NI, mat(
                ['norm_d(' + MM_1 + ', 0)', 'norm_d(' + MM_1 + ', 1)', 'norm_d(' + MM_1 + ', 2)',
                 'norm_d(' + MM_1 + ', 3)', 'norm_d(' + MM_1 + ', 4)', 'norm_d(' + MM_1 + ', 5)',
                 'norm(' + MM_1 + ', ord=2)', 'norm(' + MM_1 + ', ord=1)',
                 'norm_d(' + MM_2 + ', 0)', 'norm_d(' + MM_2 + ', 1)', 'norm_d(' + MM_2 + ', 2)',
                 'norm_d(' + MM_2 + ', 3)', 'norm_d(' + MM_2 + ', 4)', 'norm_d(' + MM_2 + ', 5)',
                 'norm(' + MM_2 + ', ord=2)', 'norm(' + MM_2 + ', ord=1)',
                 'norm_d(' + MM_3 + ', 0)', 'norm_d(' + MM_3 + ', 1)', 'norm_d(' + MM_3 + ', 2)',
                 'norm_d(' + MM_3 + ', 3)', 'norm_d(' + MM_3 + ', 4)', 'norm_d(' + MM_3 + ', 5)',
                 'norm(' + MM_3 + ', ord=2)', 'norm(' + MM_3 + ', ord=1)',
                 'norm_d(' + MM_4 + ', 0)', 'norm_d(' + MM_4 + ', 1)', 'norm_d(' + MM_4 + ', 2)',
                 'norm_d(' + MM_4 + ', 3)', 'norm_d(' + MM_4 + ', 4)', 'norm_d(' + MM_4 + ', 5)',
                 'norm(' + MM_4 + ', ord=2)', 'norm(' + MM_4 + ', ord=1)'])]
    return name_NI


def name_ni_s():
    name_MP = ['M_wei', 'M_ele', 'M_ion', 'M_rud', 'M_oue', 'M_oue/M_noe', 'M_bra']
    #    name_MP=['M_wei','M_ele','M_ion','M_rud','M_oue','M_noe','M_bra']
    name_MS = ['M_S^e', 'M_S_A^e', 'M_S_AB^e', 'M_S_ABC^e', 'np.power(M_S_bon_[M_bon_^e', 'M_S_ABC_aro_^e',
               'M_S_ABC_cyc_^e', 'M_bon__cyc_^e']
    name_MP_H = ['M_wei_H', 'M_ele_H', 'M_ion_H', 'M_rud_H', 'M_oue_H', 'M_oue_H/M_noe_H', 'M_bra_H']
    #    name_MP_H=['M_wei_H','M_ele_H','M_ion_H','M_rud_H','M_oue_H','M_noe_H','M_bra_H']
    name_MS_H = ['M_S^e', 'M_S_A^e', 'M_S_AB^e', 'M_S_ABC^e', 'M_bon_H^e', 'M_S_ABC_aro_^e', 'M_S_ABC_cyc_^e',
                 'M_bon__cyc_^e']
    name_NI = mat([])
    name_NI = name_ni(name_MS, name_MP, name_NI)
    name_NI = name_ni(name_MS_H, name_MP_H, name_NI)
    return name_NI


def ni_s(MSI):
    n_atom = len(MSI['M_atom'])
    M_S_A = MSI['M_S_A']  # 相邻矩阵
    # M_S_L=MSI['M_S_L']      #拉普拉斯矩阵
    M_S = MSI['M_S']  # 步长矩阵
    M_aro = MSI['M_aro']  # 苯环原子标1
    M_cyc = MSI['M_cyc']  # 环状原子标1
    M_cyc = M_cyc - M_aro  # 去掉苯环
    M_aro = MSI['M_aro']
    M_bon_ = MSI['M_bon_']
    M_S_B = mat(np.zeros((n_atom, n_atom)))
    M_S_B[M_S == 2] = 2  # 相间矩阵
    M_S_C = mat(np.zeros((n_atom, n_atom)))
    M_S_C[M_S == 3] = 3  # 相跳矩阵
    M_S_AB = M_S_A + M_S_B
    M_S_ABC = M_S_A + M_S_B + M_S_C
    M_S_all = {'M_S': M_S, 'M_S_A': M_S_A, 'M_S_B': M_S_B, 'M_S_C': M_S_C, 'M_S_AB': M_S_AB, 'M_S_ABC': M_S_ABC,
               'M_S_aro': np.multiply(M_S, M_aro.dot(M_aro.T)), 'M_S_ABC_aro_': np.multiply(M_S_ABC, M_aro.T),
               'M_S_cyc': np.multiply(M_S, M_cyc.dot(M_cyc.T)), 'M_S_ABC_cyc_': np.multiply(M_S_ABC, M_cyc.T),
               'M_S_A_cyc_': np.multiply(M_S_A, M_cyc.T), 'M_bon__cyc_': np.multiply(M_bon_, M_cyc.T)
               }  # 所有步长计算的矩阵,加了M_S_L的
    M_P_all = dict(m_property(MSI['M_atom']), **{'M_aro': MSI['M_aro']})  # 所有性质计算的矩阵

    in_H = np.where(M_P_all['M_wei'] == 1.00794)[0]
    M_P_H_all = h_delete(in_H, copy.deepcopy(M_P_all))  # 去 H
    M_S_r = mat_inv(copy.deepcopy({'M_S': M_S}))  # 步长矩阵取倒数
    M_S_all_f = dict(M_S_all, **{'M_S_r': copy.deepcopy(M_S_r['M_S'])})  # 合并字典
    M_S_H_all_f = h_delete(in_H, copy.deepcopy(M_S_all_f))
    M_S_bon_ = {'M_bon_': M_bon_}
    M_S_bon_H = h_delete(in_H, copy.deepcopy(M_S_bon_))

    e = np.exp(1)
    #    MP=[M_P_all['M_wei'],M_P_all['M_ele'],M_P_all['M_ion'],M_P_all['M_rud'],M_P_all['M_oue'],M_P_all['M_noe'],sum(M_S_all_f['M_S_A']).T]
    MP = [M_P_all['M_wei'], M_P_all['M_ele'], M_P_all['M_ion'], M_P_all['M_rud'], M_P_all['M_oue'],
          M_P_all['M_oue'] / M_P_all['M_noe'], sum(M_S_all_f['M_S_A']).T]  # 性质矩阵
    MS = [np.power(M_S_all_f['M_S'], e), np.power(M_S_all_f['M_S_A'], e), np.power(M_S_all_f['M_S_AB'], e),
          np.power(M_S_all_f['M_S_ABC'], e), np.power(M_S_bon_['M_bon_'], e),
          np.power(M_S_all_f['M_S_ABC_aro_'], e), np.power(M_S_all_f['M_S_ABC_cyc_'], e),
          np.power(M_S_all_f['M_bon__cyc_'], e)]  # 步长矩阵
    #    MP_H=[M_P_H_all['M_wei'],M_P_H_all['M_ele'],M_P_H_all['M_ion'],M_P_H_all['M_rud'],M_P_H_all['M_oue'],M_P_H_all['M_noe'],sum(M_S_H_all_f['M_S_A']).T]
    MP_H = [M_P_H_all['M_wei'], M_P_H_all['M_ele'], M_P_H_all['M_ion'], M_P_H_all['M_rud'], M_P_H_all['M_oue'],
            M_P_H_all['M_oue'] / M_P_H_all['M_noe'], sum(M_S_H_all_f['M_S_A']).T]
    MS_H = [np.power(M_S_H_all_f['M_S'], e), np.power(M_S_H_all_f['M_S_A'], e), np.power(M_S_H_all_f['M_S_AB'], e),
            np.power(M_S_H_all_f['M_S_ABC'], e), np.power(M_S_bon_H['M_bon_'], e),
            np.power(M_S_H_all_f['M_S_ABC_aro_'], e), np.power(M_S_H_all_f['M_S_ABC_cyc_'], e),
            np.power(M_S_H_all_f['M_bon__cyc_'], e)]
    NI = mat([])
    NI = cal_ni(MS, MP, NI)
    NI = cal_ni(MS_H, MP_H, NI)
    nnH = len(in_H)
    mw = np.sum(abs(M_P_all['M_wei']), axis=0)[0, 0]
    st_m = np.max(M_S_H_all_f['M_S'])
    st_s = np.sum(M_S_H_all_f['M_S'])
    ms = [n_atom, n_atom - nnH, mw, st_m, st_s]
    return NI, ms


def cal_Ps(gjf, Tb, T):
    try:
        warnings.filterwarnings("ignore")
        # orimsi_o = open(gjf, 'r', encoding='utf8')
        # orimsi_all = orimsi_o.readlines()
        orimsi_all=gjf
        MSI_gjf = msi_gjf(orimsi_all)
        # 关闭文件
        # orimsi_o.close()
        NIS, ms = ni_s(MSI_gjf)
        X_S = mat(NIS)
        X_M_S = mat(ms)
        ii_S = np.load('./myapp/utils/SaturatedVaporPressure/ii_left_S.npy')
        X_S = X_S[:, ii_S]
        XT_S_d = np.c_[
            copy.deepcopy(X_S), copy.deepcopy(X_S) / X_M_S[:, 1], copy.deepcopy(X_S) / np.sqrt(X_M_S[:, 2]), copy.deepcopy(
                X_S) / np.sqrt(X_M_S[:, 3]), copy.deepcopy(X_S) / np.sqrt(X_M_S[:, 4]), copy.deepcopy(X_S) / np.sqrt(
                X_M_S[:, 3] / X_M_S[:, 1])]
        re_all = np.load('./myapp/utils/SaturatedVaporPressure/result_ad_AB.npz', allow_pickle=True)
        ii_all = re_all['ii_all']
        bb_all = re_all['bb_all']
        cy_call = 40
        A_cul = bb_all[cy_call - 1][0][0]
        for i in range(0, cy_call):
            m = XT_S_d[:, ii_all[cy_call - 1][i]] * bb_all[cy_call - 1][0][i + 1]
            A_cul = A_cul + m
        P_cul = np.multiply(A_cul, (1 / (Tb - 20) - 1 / (T - 20)))
        P_cul = 101.325 * np.exp(P_cul)
        return P_cul
    except Exception as e:
        return str(e)


if __name__ == '__main__':
    cal_value = cal_Ps('./myapp/utils/SaturatedVaporPressure/504_60_9.gjf', 315, 300)
    print(cal_value)

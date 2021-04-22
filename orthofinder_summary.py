import matplotlib.pyplot as plt
import numpy as np

# 获得一张基因家族保守程度的柱状图
def get_a_normal_picture():
    with open('Orthogroups.GeneCount.csv', 'r') as f:
        ortho2num = {}
        for line in f.readlines()[1:]:
            line = line[:-1].split()
            ortho2num[line[0]] = sum([i != '0' for i in line[1:]])
        ortho2num = sorted(ortho2num.items(), key=lambda item: item[1], reverse=True)[:200]

        plt.rcParams['font.family'] = ['SimHei']              # 解决不能输出中文的问题。不区分大小写，即SimHei’效果等价于‘simhei’，中括号可以不要
        plt.rcParams['figure.autolayout'] = True              # 解决不能完整显示的问题（比如因为饼图太大，显示窗口太小）
        plt.bar(np.arange(len(ortho2num)), [item[1]-ortho2num[-1][1]+1 for item in ortho2num], color='SkyBlue')
        plt.xticks(fontsize=12, rotation=80) # 这个是改x轴标签的角度的
        plt.yticks(np.arange(1,ortho2num[0][1]-ortho2num[-1][1]+2),labels=range(ortho2num[-1][1],ortho2num[0][1]+1)) # 改刻度和刻度对应的标签
        plt.title('the largest orthogroup')
        plt.ylabel('number of species')
        plt.xlabel('orthogroup')
        plt.savefig('orthogroup_to_species.png',bbox_inches='tight')
    return

# 获得一个txt文件，里面包括特异基因数量，保守基因数量，两个以上保守基因数量，物种数量，属的数量
def the_main_summary():
    with open('orthofinder_summary_out.txt', 'w') as w1,\
        open('conserved_orthogroups.txt', 'w') as w2,\
        open('Orthogroups_UnassignedGenes.csv', 'r') as f1,\
        open('Orthogroups.GeneCount.csv', 'r') as f2:

        f2_lines = f2.readlines()
        w2.write(f2_lines[0])
        temp_num = 0
        for line in f2_lines[1:]:
            line_list = line[:-1].split('\t')
            if not '0' in line_list:
                temp_num += 1
                w2.write(line)

        w1.write('Number of specific orthogroup: {}\n'.format(len(f1.readlines())-1))
        w1.write('Number of orthogroup conserved in two species: {}\n'.format(len(f2_lines)-1))
        w1.write('Number of orthogroup conserved in all species: {}\n'.format(temp_num))
    return

if __name__=='__main__':
    get_a_normal_picture()
    the_main_summary()


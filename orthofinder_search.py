import getopt, sys


def read(file_name):
    with open(file_name, 'r') as f:
        return f.read().split()

class search():
    # 搜索物种或者同源家族
    # ------------------------------------------------------------------------------
    def __init__(self):
        pass

    # -----------------------------------------------------------------------------------
    # 这是搜索基因家族
    # 用来搜索csv文件
    def csv(self, target, file_name='Orthogroups.csv', nrows=9999999, skiprows=0, is_save = True):
        lis = []
        with open(file_name, 'r') as f, open('ortholog_to_gene.txt', 'w') as w:
            lis.append(f.readline())
            for line in f.readlines()[skiprows:nrows+skiprows]:
                temp_lis = line[:-1].replace(',',' ').split()
                for each in target:
                    if each in temp_lis:
                        w.write('{}\t{}\n'.format(temp_lis[0], each))
                        lis.append(line)
        if is_save:
            with open('search_out.txt', 'w') as w:
                w.writelines(lis)
        return lis
    # 搜索txt文件
    def txt(self, target, file_name='Orthogroups.txt', nrows=9999999, skiprows=0, is_save = True):
        lis = []
        with open(file_name, 'r') as f, open('ortholog_to_gene.txt') as w:
            for line in f.readlines()[skiprows:nrows+skiprows]:
                og = line[:-1].split(':')
                og_seq = og[1].split()
                for each in target:
                    if each in og_seq:
                        w.write('{}\t{}\n'.format(og[0], each))
                        lis.append(line)
        if is_save:
            with open('search_out.txt', 'w') as w:
                w.writelines(lis)
        return lis
    # 主函数
    def ortholog(self, target=['ESA42661'], file_name='Orthogroups.txt', nrows=9999999, skiprows=0, is_save = True):
        # 获取扩展名，来判断使用那个函数，因此该文件的需要包含'.'开头的扩展名
        suffix = file_name.rsplit('.')[-1] 
        # --------------------------------------------------------------
        if suffix == 'txt':
            return self.txt(target, file_name, nrows, skiprows, is_save)
        elif suffix == 'csv':
            return self.csv(target, file_name, nrows, skiprows, is_save)
        else:
            print('ERROR: the suffix of file name should be \'.txt\' or \'.csv\'')
            return None

    # ---------------------------------------------------------------------------------------
    # 这是对物种的筛选和搜索
    # 使用SpeciesIDs.txt获得与target对应的正确物种名
    def check_name(self, target):
        result = []
        with open('SpeciesIDs.txt', 'r') as f:
            lines = f.readlines()
            for line in lines:
                for each in target:
                    if line.find(each) != -1:
                        result.append(line.rsplit('.fa')[0].split()[1])
        # print(result)
        return result
    # 获得索引
    def search_index(self, target, species):
        species_dict = {} # 用字典来加速搜索
        for index in range(len(species)):
            species_dict[species[index]] = index
        # print(species_dict)
        try:
            return [species_dict[each] for each in target]
        except KeyError as error:
            print('Error: species name not right')
    # 主函数
    def species(self, target=[], file_name='search_out.txt', file_list = None):
        # 检查物种名是否正确
        try:
            right_target = self.check_name(target)
        except:
            right_target = target
            print('Wrong: species may not right')
        # ------------------
        if file_list == None:
            with open(file_name, 'r') as f:
                lines = f.readlines()
        else:
            lines = file_list
            # 写入首行
        # ------------------
        with open('run_out.tsv', 'w') as w:
            w.write('\t{}\n'.format('\t'.join(right_target)))
            # 得到每个物种的索引index
            species = lines[0]
            species = species.split('\t')
            indexs = self.search_index(right_target, species)
            # print(right_target, species)
            # 在基因家族中筛选
            lines = lines[1:]
            for line in lines:
                line = line.split('\t')
                # print(line)
                write_line = [line[0]] # 用来临时存储写入的一行
                for index in indexs:
                    write_line.append(line[index])
                w.write('\t'.join(write_line))
                w.write('\n')

    # ----------------------------------------------------------------------------------
    # 物种和同源组同时搜索
    def all(self, ortholog = [], species = [], file_name = 'Orthogroups.csv'):
        ortholog_result = self.ortholog(ortholog, is_save=False, file_name=file_name)
        # print(ortholog_result)
        return self.species(target = species, file_list = ortholog_result)
        

if __name__ == '__main__':
    # 获得参数
    opts, args = getopt.getopt(sys.argv[1:], 'g:s:f:')
    geneid_filename = None
    species_filename = None
    search_filename = 'Orthogroups.csv'
    for o, a in opts:
        if o in ('-g'):
            geneid_filename = a
        if o in ('-s'):
            species_filename = a
        if o in ('-f'):
            search_filename = a
    # ---------------------
    s =search()
    if geneid_filename == None and species_filename == None:
        print('Error: no parameter')
    elif geneid_filename != None and species_filename == None:
        s.ortholog(read(geneid_filename), search_filename)
    elif geneid_filename == None and species_filename != None:
        s.species(read(species_filename), search_filename)
    else:
        s.all(read(geneid_filename), read(species_filename), search_filename)
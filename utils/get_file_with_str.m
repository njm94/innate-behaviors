function f = get_file_with_str(data_dir, str_in)
file_list = getAllFiles(data_dir);
file_idx = contains(file_list, str_in);
f = file_list{file_idx};
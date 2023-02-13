"Check if `dirname` exists, if not create the corresponding directory."
function check_dir(dirname)
    isdir(dirname) || mkdir(dirname)
end

#!python
""" metafil.py

    Dida Markovic, 29/03/2017

    Collection of old filenaming functions and other dinky meta tools.
    Many are from forumns etc (links in comments).
"""
import inspect
import glob
import time
import subprocess
import os
import os.path
import unicodedata
import pkg_resources

""" ---------------------- DEBUGGING TOOLS ----------------------- """


def lineno(string=''):
    # Returns the current line number in our program. How come not always the line in this file??
    if string != '':
        return "L." + str(inspect.currentframe().f_back.f_lineno) + ': ' + str(string)
    else:
        return "L." + str(inspect.currentframe().f_back.f_lineno)
# http://code.activestate.com/recipes/145297-grabbing-the-current-line-number-easily/


""" ------------------- END OF DEBUGGING TOOLS ------------------- """

""" ---------------------- FILENAME TOOLS ------------------------ """


def ts(t=(2015, 3, 1, 0, 0, 0, 0, 0, 0)):
    return str(int(time.time() - time.mktime(t)))


def validate_filename(filepath):
    """ Raises an Exception if the file does not exist.

        Parameters
        ----------
        filepath : str
            input filename

        Returns
        -------
        filepath : str
            validated, absolute path to the file
    """
    if not os.path.isfile(filepath):
        raise Exception('File not found: ' + filepath)
    return os.path.abspath(filepath)


def new_suffix(filename, suffix):
    """ If fileout exists in the file system it adds [...].v#, where # is the number that
           creates a filename that doesn't exist yet.

    Parameters
    ----------
    filename  : str
            input filename (can be with path) with undesired suffix (defined as after last dot!)
    suffix  : str
            new suffix (dot separator added automatically, so don't include here!)

    Returns
    -------
    filename : str
        filename with new suffix
    """

    # Read the suffix, which is assumed to be a string after the last dot
    return '.'.join(filename.split('.')[:-1] + [suffix])


def increment_filename(fileout):
    """ If fileout exists in the file system it adds [...].v#, where # is the number that
           creates a filename that doesn't exist yet.

    Parameters
    ----------
    fileout  : str
            full desired pathname (will be incremented if exists)

    Returns
    -------
    filepath : str
        new path not associated with a existing file
    """

    # Read the suffix, which is assumed to be a string after the last dot
    suffix = '.' + fileout.split('.')[-1]

    # If it already exists, increment it
    n = 1
    while os.path.exists(fileout):
        n += 1
        if n == 2:
            fileout = '.'.join(fileout.split('.')[:-1]) + '.v' + str(n) + suffix
            continue
        fileout = '.'.join(fileout.split('.')[:-2]) + '.v' + str(n) + suffix

    # Return the absolute path so there are no relative path problems later
    return os.path.abspath(fileout)


def _ensurelist(supposed_list):
    if isinstance(supposed_list, list):
        return supposed_list
    else:
        return [supposed_list]


def file_list(inpath):
    """ Returns a list of strings that are valid paths in the system.

    Parameters
    ----------
    inpath : str or list
        string path or a list of them, can be either a single file,
        single folder or a list of files or folders

    Returns
    -------
    file_list : list
        a list of strings, each is a valid path to a file in the system
    """
    file_list = []
    for ip in _ensurelist(inpath):
        if os.path.isfile(ip):
            file_list.append(ip)
        elif os.path.isdir(ip):
            file_list.extend(glob.glob(ip + '/*'))
    return file_list


def strdiff(filenames, spliton='_'):
    """ This finds the part of all the strings that differs between
        the strings. It may take long if have a lot of strings, but I am not sure.
        The problem here is that order is not retained well... Annoying.
        Note also that if you have the the string that is unique to each filename
        also appearing in the path, you need to run securediff on the result.
        Finally, note that this doesn't work if many parts of the filename are varying
        if they are separated by _, . or -! """

    if len(filenames) == 0:
        raise Exception("There are no files in the input to diff!")

    thingstospliton = ['_', '.', '-']

    # Stop if last separator in list
    try:
        new_spliton = thingstospliton[thingstospliton.index(spliton) + 1]
    except IndexError:
        new_spliton = ''

    # Make a dictionary of sets of substrings with full filenames as keys
    namelist = {}
    for filename in filenames:
        namelist[filename] = filename.split('/')[-1].split(spliton)

    # Find intersections of all the other sets with the first, keep common strings
    for filename in filenames:
        if filename is filenames[0]:
            saveset = namelist[filename]
            continue
        saveset = list(set(saveset) & set(namelist[filename]))

    # Remove common substrings and re-merge sets to strings
    for filename in filenames:
        namelist[filename] = new_spliton.join(set(namelist[filename]) - set(saveset))

    # Nest this function - try to split it on '.' and '-' as well
    if new_spliton != '':
        tmp = strdiff(namelist.values(), spliton=new_spliton)
        # Populate the dictionary to return
        for filename in filenames:
            namelist[filename] = tmp[namelist[filename]]

    return namelist


def securediff(fdict):
    """ This adds some security padding around the unique part of a string to make it less
        likely that there will be repetitions of it in the string. E.g. if it's just a 0,
        and have numbers in the path too.
        Note that there should be only 1 unique string and up to 2 common strings
        surrounding it!"""

    if len(fdict) == 0:
        raise Exception("There are no files in the input to diff!")

    common = None
    problematic = []
    for fname, uniqstr in fdict.items():

        # Find the string in common
        if fname.count(uniqstr) == 1:
            fparts = fname.split(uniqstr)

            # Check that the common string is the same in all the filenames
            if common is not None:
                for string in common:
                    if string not in fparts:
                        raise Exception('Something has gone terribly wrong!')
            common = fparts

        # List the files where the characteristic string repeats in the path
        elif fname.count(uniqstr) > 1:
            problematic.append(fname)

        # Check that it does appear in the path
        elif fname.count(uniqstr) == 0:
            print(fname, uniqstr)
            raise Exception('Something has gone terribly wrong!')

    # Check that common string found is also present in the problematic filenames
    for problem in problematic:
        for string in common:
            if string not in problem:
                raise Exception('Something has gone terribly wrong!')

    # Assuming unique string is in the middle of two (possibly empty) common
    #   strings, pad the unique string so that it always appears only once.
    if len(problematic) > 0:

        for fname, uniqstr in fdict.items():

            if len(common[0]) > 0:
                new_uniqstr = common[0][-1] + uniqstr
            if len(common[1]) > 0:
                new_uniqstr += common[1][0]

            fdict[fname] = new_uniqstr

        fdict = securediff(fdict)

    return fdict


def fnamediff(pathpattern, padding=''):
    """ This finds the part of all the filenames that differs between
        the files. It may take long if have a lot of files, but I am not sure."""

    # this is complete non-sense:
    if '*' not in pathpattern:
        if pathpattern[-1] != '/':
            pathpattern += '/'
        pathpattern += '*'

    # if len([pathpattern.start() for m in re.finditer('*', pathpattern)])>1:
    #   raise TooSpecificErr()

    splitname = pathpattern.split('*')

    prefix = splitname[0]
    lenpre = len(prefix)

    if len(splitname) > 1:
        suffix = splitname[-1]
    else:
        suffix = ''
    lensu = len(suffix)

    filenames = glob.glob(pathpattern)
    nfiles = len(filenames)

    refname = filenames[0]

    midtix = refname.replace(prefix, '').replace(suffix, '')
    # Prefix
    for letter in midtix:
        new_filenames = glob.glob(prefix + letter + '*' + suffix)
        if len(new_filenames) == nfiles:
            prefix = prefix + letter

    # Suffix
    for letter in midtix[::-1]:
        new_filenames = glob.glob(prefix + '*' + letter + suffix)
        if len(new_filenames) == nfiles:
            suffix = letter + suffix

    i = 0
    names = {}
    for filename in filenames:
        names[filename] = padding + filename.replace(prefix, '').replace(suffix, '')
        i += 1

    return names


""" ------------------- END OF FILENAME TOOLS -------------------- """

""" ------------------------- GIT TOOLS -------------------------- """


def searchup(path, filename, maxstep=3):
    """ Recursively searches up the path until it finds filename.
        Traverses max maxstep levels.
    """
    if maxstep < 0:
        raise IOError('Path not found!')
    path = os.path.abspath(path)
    full_filename = os.path.join(path, filename)
    if os.path.exists(full_filename):
        return full_filename
    else:
        return searchup(path + '/..', filename, maxstep - 1)

# A class that contains the Git environment at the time of it's initialisation.
# Currently it uses the subprocess module to speak to Git through the system.
# Ideally some day it would use the GitPython module or Dulwich or something like it.


class GitEnv(object):

    def __init__(self, home='.', name=None):
        home = os.path.dirname(os.path.realpath(home))
        self.name = name
        try:
            self.git_dir = searchup(home, '.git')
        except IOError:
            self.git_dir = home
            self.isrepo = False
        else:
            self.isrepo = True
        self.hash, self.author, self.date = [str(s) for s in self.get_commit()]
        self.url = str(self.get_remote()).split('@')[-1]
        self.branch = str(self.get_branch())
        self.repo = str(self.get_repo())
        self.printstart = ''
        self.version = str(self.get_version(name))

    def __str__(self):
        startline = self.printstart
        as_string = startline + "This was generated at " + time.strftime('%H:%M:%S, %d %b %Y')
        as_string += "\n" + startline + "\t by"
        if self.name is not None:
            as_string += " the " + self.name + " code, version " + self.version + ", from"
        else:
            as_string += " code from"
        if self.isrepo:
            as_string += " the Git repo:"
            as_string += "\n" + startline + "\t\t" + self.url + ","
            as_string += "\n" + startline + "\t on the " + self.branch + " branch,"
            as_string += "\n" + startline + "\t with commit: " + self.hash[:10]
            if self.is_dirty():
                as_string += "-dirty"
            as_string += "\n" + startline + "\t\t from " + self.date + ", "
            as_string += "\n" + startline + "\t\t by " + self.author
        else:
            as_string += " the local folder (no Git repo found):"
            as_string += "\n" + startline + "\t\t" + self.url + ","
            as_string += "\n" + startline + "\t by " + self.author
        as_string += "."
        return unicodedata.normalize('NFKC', as_string.decode("unicode-escape")).encode('ascii', 'ignore')

    def set_print(self, startline):
        self.printstart = startline

    def get_git_cmd(self, args=[]):
        cmd = ['git']
        if self.git_dir is not None:
            cmd.append('--git-dir')
            cmd.append(self.git_dir)
        for one in args:
            cmd.append(one)
        return cmd

    def get_hash(self, nochar=7, sep=''):
        return sep + self.hash[0:nochar] + sep

    def describe(self):
        if not self.isrepo:
            raise ValueError("Not a git repo: %s" % self.git_dir)
        cmd = subprocess.Popen(self.get_git_cmd(['describe', ]), stdout=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        return cmd_out.strip()

    # Get the hash, author and date of the most recent commit of the current repo.
    def get_commit(self):
        if not self.isrepo:
            return [time.gmtime(), os.getenv('USER'), time.strftime('%d/%b/%Y', time.gmtime())]
        cmd = subprocess.Popen(self.get_git_cmd(['log', '-n', '1']), stdout=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        newlist = []
        for entry in cmd_out.strip().split('\n'):
            if entry == '':
                continue
            entry = entry.split(' ')
            # This is a hack, should use a dict so can be sure what we are reading in:
            if 'commit' in entry[0] or 'Author' in entry[0] or 'Date' in entry[0]:
                newlist.append(' '.join(entry[1:]).strip())
        if len(newlist) != 3:
            raise Exception('No commit found!')
        return newlist

    # At the moment this only gives the first url in what git returns.
    # Eventually it'd be nice if you could get the origin url, the fetch...
    def get_remote(self):
        if not self.isrepo:
            return os.path.abspath(self.git_dir)
        cmd = subprocess.Popen(self.get_git_cmd(['remote', '-v']), stdout=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        if bool(cmd_out):
            try:
                return cmd_out.strip().split('https://')[1].split(' ')[0]
            except IndexError:
                ssh_url = cmd_out.strip().split('git@')[1].split(' ')[0]
                return ssh_url.replace(':', '/')
        else:
            return None

    def get_branch(self):
        if not self.isrepo:
            return None
        cmd = subprocess.Popen(self.get_git_cmd(['branch']), stdout=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        branches = cmd_out.strip().splitlines()
        for branch in branches:
            if '*' in branch:
                return branch.replace('*', '').strip()

    def get_repo(self):
        if not self.isrepo and self.name is None:
            return os.path.basename(self.git_dir)
        elif not self.isrepo:
            return self.name
        cmd = subprocess.Popen(self.get_git_cmd(['rev-parse', '--show-toplevel']), stdout=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        repo = cmd_out.strip().split('/')[-1]
        if self.name is not None:
            assert self.name == repo, \
                "Misatch between passed distribution name (" + str(self.name) + ") and repo name (" + str(repo) + ")!"
        return repo

    def is_dirty(self):
        if not self.isrepo:
            return None
        cmd = subprocess.Popen(self.get_git_cmd(['status']), stdout=subprocess.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        return 'modified' in cmd_out

    def get_version(self, name=None):
        if name is None:
            return None
        else:
            try:
                return pkg_resources.require(name)[0].version
            except pkg_resources.DistributionNotFound:
                raise ImportError("A '" + name +
                                  "' distribution was not found, but installation is required to access its version.")


""" ---------------------- END OF GIT TOOLS ---------------------- """

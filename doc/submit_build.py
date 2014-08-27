import subprocess

# Get SHA1 of current commit (after pulling)
p = subprocess.Popen(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE)
so = p.communicate()
current_hash = so[0].strip()

# Get SHA1 of last-built commit
with open('last_docs_built.txt', 'r+') as f:
    last_hash = f.readline()

    if current_hash == last_hash:
        print 'Docs are up to date.'
    else:
        print 'Building docs...'
        bsub_args = ['bsub', '-q', 'short', '-W', '6:00', 'make', 'clean',
                     'html']
        subprocess.call(bsub_args)
        f.seek(0)
        f.write(current_hash)
        f.truncate()

import pdfkit, subprocess

'''
urls=[]
g = open('url.txt','w+')
with open('DataStructures.html','r+') as f:
	urls=map(lambda x:"http://en.wikipedia.org"+x[13:x.index('" title')],filter(lambda x:x.startswith('<li><a href="/wiki'),f.readlines()))
g.write(str(urls))
g.close()
'''

with open('url.txt','r') as g:
    c =  g.readlines()[0].split()
    for a in c:
        b = a.strip('[]\',')
        c = b[b.index('wiki/')+5:]
        subprocess.call(["wkhtmltopdf",b,c+'.pdf'])

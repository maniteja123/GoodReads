import pdfkit, subprocess

'''
urls=[]
g = open('url.txt','w+')
with open('DataStructures.html','r+') as f:
	urls=map(lambda x:"http://en.wikipedia.org"+x[13:x.index('" title')],filter(lambda x:x.startswith('<li><a href="/wiki'),f.readlines()))
g.write(str(urls))
g.close()
'''

c = 'http://en.wikipedia.org/wiki/(a,b)-tree'
with open('url.txt','r') as g:
    for a in g.readlines()[0].split():
        b = a.strip('[]')
        c = b[b.index('wiki/')+5:-2]
        #pdfkit.from_url(b, c+'.pdf')
        subprocess.call(["wkhtmltopdf",b,c+'.pdf'])

    


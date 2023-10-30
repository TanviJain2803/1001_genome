#in order to edit github files use the following steps to authenticate your ssh

#log in to ssh
ssh-keygen -t ecdsa -b 521 -C "tanvi.jain2803@gmail.com"
cat /gale/netapp/home/jswift/.ssh/id_ecdsa.pub
#add key to Github
ssh -T git@github.com
eval `ssh-agent -s`
ssh-add ~/.ssh/id_ecdsa
git push

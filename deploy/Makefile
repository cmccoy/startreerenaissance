PK := $(HOME)/.ssh/id_rsa
KEY_ID := cmccoy@stoat
REGION := us-west-2
ZONE := $(REGION)c
INSTANCE_TYPE := c3.2xlarge
MASTER_INSTANCE_TYPE := m1.large
SPOT_BID := 0.61
CLUSTER_NAME := test-spark
SLAVES ?= 6

ANSIBLE_NOCOWS=1
export ANSIBLE_NOCOWS

get-master:
	spark-ec2/spark-ec2 --region $(REGION) get-master $(CLUSTER_NAME)

ssh:
	spark-ec2/spark-ec2 --identity-file $(PK) --region $(REGION) login $(CLUSTER_NAME)

launch:
	spark-ec2/spark-ec2 \
		--identity-file $(PK) \
		--key-pair $(KEY_ID) \
		--slaves $(SLAVES) \
		--spot-price $(SPOT_BID) \
		--region $(REGION) \
		--zone $(ZONE) \
		--instance-type $(INSTANCE_TYPE) \
		--master-instance-type $(MASTER_INSTANCE_TYPE) \
		--wait 240 \
		launch $(CLUSTER_NAME) && \
		python tag-instances.py $(CLUSTER_NAME)

destroy:
	spark-ec2/spark-ec2 \
		--region $(REGION) \
		destroy $(CLUSTER_NAME)

login:
	spark-ec2/spark-ec2 \
		-i $(PK) \
		-k $(KEY_ID) \
		--region $(REGION) \
		login $(CLUSTER_NAME)

ansible-list:
	ansible-playbook -i ec2.py site.yml --list-hosts

ansible-provision:
	ANSIBLE_NOCOWS=1 ansible-playbook -i ec2.py site.yml

ansible-clear-old-work:
	ansible -i ec2.py -u ec2-user --sudo -m shell -a 'ls -d /root/spark/work/app-* | head -n -1 | xargs rm -r' tag_spark_node_type_slave

.PHONY: launch get-master login ansible-list ansible-provision ansible-clear-old-work destroy ssh

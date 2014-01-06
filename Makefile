PK := $(HOME)/.ssh/id_rsa
KEY_ID := cmccoy@stoat
REGION := us-west-2
INSTANCE_TYPE := m1.large
CLUSTER_NAME := test-spark
SLAVES ?= 1

get-master:
	spark-ec2/spark-ec2 --region $(REGION) get-master $(CLUSTER_NAME)

launch:
	spark-ec2/spark-ec2 \
		--identity-file $(PK) \
		--key-pair $(KEY_ID) \
		--slaves $(SLAVES) \
		--spot-price 0.09 \
		--region $(REGION) \
		--instance-type $(INSTANCE_TYPE) \
		--wait 300 \
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
	ansible-playbook -i ec2.py site.yml

.PHONY: launch get-master login ansible-list ansible-provision destroy

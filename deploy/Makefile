PK := $(HOME)/.ssh/id_rsa
KEY_ID := cmccoy@stoat
REGION := us-west-2
ZONE := $(REGION)b
INSTANCE_TYPE := c3.8xlarge
MASTER_INSTANCE_TYPE := m3.large
SPOT_BID ?= 1.26
CLUSTER_NAME ?= test-spark
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
		python tag-instances.py $(CLUSTER_NAME) && \
		python add-cloudwatch-alarm.py $(CLUSTER_NAME)

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
	ANSIBLE_NOCOWS=1 ansible-playbook -f 9 -i ec2.py site.yml --limit tag_spark_cluster_name_$(CLUSTER_NAME)

ansible-clear-old-work:
	ansible -i ec2.py -u ec2-user --sudo -m shell -a 'ls -d /root/spark/work/app-* | head -n -1 | xargs rm -r' tag_spark_node_type_slave

ansible-restart-slaves:
	ansible -i ec2.py -u ec2-user --sudo -m shell -a '/root/spark/bin/stop-slaves.sh; /root/spark/bin/start-slaves.sh' tag_spark_node_type_master

view-history:
	aws ec2 describe-spot-price-history --max-items 100 --instance-types $(INSTANCE_TYPE) --product-descriptions Linux/UNIX --availability-zone $(ZONE)

.PHONY: launch get-master login ansible-list ansible-provision ansible-clear-old-work destroy ssh view-history



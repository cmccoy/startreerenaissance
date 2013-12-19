PK := $(HOME)/.ssh/id_rsa
KEY_ID := cmccoy@stoat
REGION := us-west-2
INSTANCE_TYPE := m1.medium
CLUSTER_NAME := test-spark

launch:
	spark-ec2/spark-ec2 \
		-i $(PK) \
		-k $(KEY_ID) \
		--spot-price 0.09 \
		--region $(REGION) \
		--instance-type $(INSTANCE_TYPE) \
		--wait 300 \
		launch $(CLUSTER_NAME)

get-master:
	spark-ec2/spark-ec2 --region $(REGION) get-master $(CLUSTER_NAME)

login:
	spark-ec2/spark-ec2 \
		-i $(PK) \
		-k $(KEY_ID) \
		--region $(REGION) \
		login $(CLUSTER_NAME)

.PHONY: launch get-master login

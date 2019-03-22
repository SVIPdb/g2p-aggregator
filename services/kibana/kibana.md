## Notes about Kibana

G2P v0.12 is using Kibana verison 6.4.3

A known issue of this version of kibana and kibana version 5.6.0, is that single count visualizations, like that used top left of the g2p dashboard, have unexpected scrolling behavior.

This issue has been [documented](https://github.com/elastic/kibana/issues/14066) to the Elastic company, and is likely to be resolved in newer versions of kibana; perhaps v5.7.0.

We intend to update to the newer version once this issue is resolved by Elastic.

import unittest

from backend.app import app


class ApiSmokeTest(unittest.TestCase):
    def setUp(self):
        self.client = app.test_client()

    def test_health(self):
        response = self.client.get("/api/health")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json(), {"status": "ok"})

    def test_inchi_layers(self):
        response = self.client.get("/api/inchi_layers")

        self.assertEqual(response.status_code, 200)
        self.assertGreater(len(response.get_json()), 0)

    def test_compare_inchis(self):
        response = self.client.post(
            "/api/compare_inchis",
            json={"inchi1": "CCO", "inchi2": "CCO"},
        )

        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.get_json()["results"]["COMPLETE_IDENTITY"])

    def test_compare_inchis_custom_accepts_levels_alias(self):
        response = self.client.post(
            "/api/compare_inchis_custom",
            json={"inchi1": "CCO", "inchi2": "CCO", "levels": ["complete_identity"]},
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(list(response.get_json()["results"].keys()), ["COMPLETE_IDENTITY"])

    def test_generate_3d_requires_input(self):
        response = self.client.post("/api/generate_3d", json={})

        self.assertEqual(response.status_code, 400)


if __name__ == "__main__":
    unittest.main()
